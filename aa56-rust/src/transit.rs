//! Cálculo de Nascer (Rise), Passagem Meridiana (Transit) e Pôr (Set) dos astros (`trnsit.c`),
//! incluindo busca por bissecção para eventos próximos a regiões polares (`no_rise_set` / `search_halve`),
//! além da tabela dedicada de nascer e pôr da Lua (`moonrise.c`).

use crate::earth_orientation::{refrac, sidrlt};
use crate::math_time::{deltat, jtocal_parts, modtp, DTR, J2000, RTD};
use std::f64::consts::PI;

use serde::Serialize;

const DISFAC: f64 = 2.3454780e4; // Raios terrestres por AU (trnsit.c)
const COSSUN: f64 = -0.014543897651582657; // cos(90° 50') para o Sol
const COSZEN: f64 = -9.8900378587411476e-3; // cos(90° 34') para estrelas e planetas

pub const HORIZON_SUN: i32 = 0;
pub const HORIZON_PLANET_STAR: i32 = 1;
pub const HORIZON_MOON: i32 = 3;

pub trait IntoRaDecDist {
    fn into_ra_dec_dist(self) -> (f64, f64, f64);
}

impl IntoRaDecDist for (f64, f64, f64) {
    fn into_ra_dec_dist(self) -> (f64, f64, f64) {
        self
    }
}

impl IntoRaDecDist for (f64, f64, f64, f64) {
    fn into_ra_dec_dist(self) -> (f64, f64, f64) {
        (self.0, self.1, 1.0)
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct RiseTransitSet {
    pub has_transit: bool,
    pub has_rise: bool,
    pub has_set: bool,
    pub jd_rise_ut: Option<f64>,
    pub jd_transit_ut: Option<f64>,
    pub jd_set_ut: Option<f64>,
    pub rise_hms_ut: String,
    pub transit_hms_ut: String,
    pub set_hms_ut: String,
    pub rise_cal_ut: String,
    pub transit_cal_ut: String,
    pub set_cal_ut: String,
    pub transit_elevation_deg: f64,
    pub visible_hours: f64,
}

pub type RiseSetTransit = RiseTransitSet;

/// Formata a fração de dia UT de um Dia Juliano `jd` (onde `.5` é 00:00 UT) como `HH:MM:SS UT`.
pub fn jd_to_hms_ut(jd: f64) -> String {
    let frac = (jd + 0.5) - (jd + 0.5).floor();
    let mut total_sec = (frac * 86400.0).round() as i64;
    if total_sec >= 86400 {
        total_sec -= 86400;
    }
    let h = total_sec / 3600;
    let m = (total_sec % 3600) / 60;
    let s = total_sec % 60;
    format!("{:02}:{:02}:{:02} UT", h, m, s)
}

/// Aproximação analítica de um passo de `trnsit.c` (linhas 34-138).
/// Retorna `(f_trnsit, r_trnsit, r_rise, r_set, elevation_threshold_deg)`.
pub fn trnsit_step(
    jd_ut: f64,
    lha_rad: f64,
    dec_rad: f64,
    dist_au: f64,
    dradt_rad_day: f64,
    ddecdt_rad_day: f64,
    glat_deg: f64,
    obj_kind: i32, // 0 = Sun, 3 = Moon, default = Planet/Star
) -> (bool, f64, f64, f64, f64) {
    let tpi = 2.0 * PI;
    let x_lat = glat_deg * DTR;
    let coslat = x_lat.cos();
    let sinlat = x_lat.sin();
    let cosdec = dec_rad.cos();
    let sindec = dec_rad.sin();

    let date0 = (jd_ut - 0.5).floor() + 0.5;
    let x = (jd_ut - date0) * tpi;

    let mut y = lha_rad;
    while y < -PI {
        y += tpi;
    }
    while y > PI {
        y -= tpi;
    }
    let lhay = y;
    let denom = -dradt_rad_day / tpi + 1.00273790934;
    let r_trnsit = x - y / denom;

    if coslat.abs() < 1e-12 || cosdec.abs() < 1e-12 {
        return (false, r_trnsit, -10.0, -10.0, -0.5667);
    }

    let (n_hor, elev_thresh) = match obj_kind {
        0 => (COSSUN, -0.8333333333333333),
        3 => {
            let n_px = (1.0 / (DISFAC * dist_au.max(1e-6))).clamp(-1.0, 1.0);
            let d_px = n_px.asin();
            let sdiam_rad = 0.2725076 * d_px + 3.874e-7;
            let n_val = (-9.890199094634534e-3 - sdiam_rad + d_px).sin();
            let thresh = -34.0 / 60.0 - sdiam_rad * RTD;
            (n_val, thresh)
        }
        _ => (COSZEN, -0.5666666666666666),
    };

    let y_cos = (n_hor - sinlat * sindec) / (coslat * cosdec);
    if y_cos <= -1.0 || y_cos >= 1.0 {
        return (false, r_trnsit, -10.0, -10.0, elev_thresh);
    }

    let mut z = -ddecdt_rad_day * (sinlat + COSZEN * sindec);
    z /= tpi * coslat * cosdec * cosdec;
    z /= (1.0 - y_cos * y_cos).sqrt();
    let y_arc = y_cos.acos();

    let r_rise = x - (lhay + y_arc) * (1.0 + z) / denom;
    let r_set = x - (lhay - y_arc) * (1.0 - z) / denom;

    (true, r_trnsit, r_rise, r_set, elev_thresh)
}

/// Busca por bissecção de intervalo (`search_halve` em `trnsit.c` linhas 461-495).
fn search_halve<F>(t1_in: f64, _y1_in: f64, t2_in: f64, y2_in: f64, elev_thresh: f64, eval_elev: &mut F) -> f64
where
    F: FnMut(f64) -> f64,
{
    let mut t1 = t1_in;
    let mut t2 = t2_in;
    let mut e2 = y2_in - elev_thresh;
    let mut tm = 0.5 * (t1 + t2);

    while (t2 - t1).abs() > 0.00001 {
        tm = 0.5 * (t1 + t2);
        let ym = eval_elev(tm);
        let em = ym - elev_thresh;
        if em * e2 > 0.0 {
            t2 = tm;
            e2 = em;
        } else {
            t1 = tm;
        }
    }
    tm
}

/// Busca de fallback quando a aproximação inicial falha em encontrar nascer/pôr (`no_rise_set` em `trnsit.c` linhas 351-455).
fn no_rise_set<F>(
    t_trnsit0: f64,
    el_trnsit0: f64,
    elev_thresh: f64,
    eval_elev: &mut F,
) -> (Option<f64>, Option<f64>)
where
    F: FnMut(f64) -> f64,
{
    let mut t_rise_out = None;
    let mut t_set_out = None;

    // Passo para trás em direção ao trânsito inferior anterior
    let mut t1 = t_trnsit0;
    let mut t2 = t1;
    let mut el1 = el_trnsit0;
    loop {
        if (t_trnsit0 - t2) > 1.0 {
            break;
        }
        t2 = t1 - 0.25;
        let el2 = eval_elev(t2);
        if el2 < elev_thresh {
            t_rise_out = Some(search_halve(t2, el2, t1, el1, elev_thresh, eval_elev));
            break;
        }
        if el2 < el1 {
            t1 = t2;
            el1 = el2;
            continue;
        }
        break;
    }

    // Passo para frente em direção ao próximo trânsito inferior
    t1 = t_trnsit0;
    t2 = t1;
    el1 = el_trnsit0;
    loop {
        if (t2 - t_trnsit0) > 1.0 {
            break;
        }
        t2 = t1 + 0.25;
        let el2 = eval_elev(t2);
        if el2 < elev_thresh {
            t_set_out = Some(search_halve(t2, el2, t1, el1, elev_thresh, eval_elev));
            break;
        }
        if el2 < el1 {
            t1 = t2;
            el1 = el2;
            continue;
        }
        break;
    }

    (t_rise_out, t_set_out)
}

/// Solução iterativa completa de Nascer, Trânsito e Pôr (`iter_trnsit` de `trnsit.c`).
/// Aceita parâmetros completos de observação e closure `eval_radec_dist(jd_tdt) -> (ra_rad, dec_rad, ...)`
pub fn compute_rise_transit_set<F, R>(
    jd_ut_start: f64,
    obs_lon_deg: f64,
    obs_lat_deg: f64,
    _geo_obs: &(f64, f64),
    obj_kind: i32,
    temp_c: f64,
    pressure_mb: f64,
    mut eval_radec_dist: F,
) -> RiseTransitSet
where
    F: FnMut(f64) -> R,
    R: IntoRaDecDist,
{
    let tpi = 2.0 * PI;
    let mut eval_at_ut = |t_ut: f64| -> (bool, f64, f64, f64, f64, f64, f64) {
        let y_ep = 2000.0 + (t_ut - J2000) / 365.25;
        let dt_sec = deltat(y_ep);
        let tdt = t_ut + dt_sec / 86400.0;
        let (ra, dec, dist) = eval_radec_dist(tdt).into_ra_dec_dist();

        let step_d = 0.001; // igual a domoon.c (TDT - 0.001)
        let (ra_p, dec_p, _) = eval_radec_dist(tdt + step_d).into_ra_dec_dist();
        let (ra_m, dec_m, _) = eval_radec_dist(tdt - step_d).into_ra_dec_dist();
        let mut dra = ra_p - ra_m;
        if dra > PI {
            dra -= tpi;
        } else if dra < -PI {
            dra += tpi;
        }
        let dradt = dra / (2.0 * step_d);
        let ddecdt = (dec_p - dec_m) / (2.0 * step_d);

        let last_sec = sidrlt(t_ut, tdt, obs_lon_deg);
        let last_rad = last_sec * DTR / 240.0;
        let lha = modtp(last_rad - ra);

        let glat_rad = DTR * obs_lat_deg;
        let sin_alt = glat_rad.sin() * dec.sin() + glat_rad.cos() * dec.cos() * lha.cos();
        let alt_true_deg = RTD * sin_alt.clamp(-1.0, 1.0).asin();
        let alt_app_deg = alt_true_deg + refrac(alt_true_deg, pressure_mb, temp_c);

        let (ok, r_tr, r_ri, r_se, elev_thresh) =
            trnsit_step(t_ut, lha, dec, dist, dradt, ddecdt, obs_lat_deg, obj_kind);
        (ok, r_tr, r_ri, r_se, alt_app_deg, alt_true_deg, elev_thresh)
    };

    // 1. Itera o horário de Trânsito Meridiano (`iter_trnsit` linhas 192-209)
    let mut t1 = jd_ut_start;
    let mut r_rise_init = -10.0;
    let mut r_set_init = -10.0;
    let mut f_trnsit = false;
    let mut elev_trnsit = 0.0;
    let mut elev_true_trnsit = 0.0;
    let mut elev_thresh = -0.5667;
    let mut date_trnsit = (t1 - 0.5).floor() + 0.5;

    for _ in 0..11 {
        let t0 = t1;
        let date = (t0 - 0.5).floor() + 0.5;
        let (ok, r_tr, r_ri, r_se, el_app, el_true, thresh) = eval_at_ut(t0);
        f_trnsit = ok;
        r_rise_init = r_ri;
        r_set_init = r_se;
        elev_trnsit = el_app;
        elev_true_trnsit = el_true;
        elev_thresh = thresh;
        date_trnsit = date;
        t1 = date + r_tr / tpi;
        if (t1 - t0).abs() <= 0.0001 {
            break;
        }
    }
    let t_trnsit = t1;
    let (_, _, _, _, _, _, tr_cal) = jtocal_parts(t_trnsit);

    if !f_trnsit {
        if elev_true_trnsit > elev_thresh {
            let mut eval_elev = |t: f64| -> f64 { eval_at_ut(t).5 };
            let (opt_r, opt_s) = no_rise_set(t_trnsit, elev_true_trnsit, elev_thresh, &mut eval_elev);
            if opt_r.is_some() || opt_s.is_some() {
                let has_r = opt_r.is_some();
                let has_s = opt_s.is_some();
                let r_hms = opt_r.map(jd_to_hms_ut).unwrap_or_else(|| "—".to_string());
                let s_hms = opt_s.map(jd_to_hms_ut).unwrap_or_else(|| "—".to_string());
                let r_cal = opt_r.map(|j| jtocal_parts(j).6).unwrap_or_else(|| "—".to_string());
                let s_cal = opt_s.map(|j| jtocal_parts(j).6).unwrap_or_else(|| "—".to_string());
                let vis = match (opt_r, opt_s) {
                    (Some(tr), Some(ts)) if ts > tr && (ts - tr) < 1.0 => 24.0 * (ts - tr),
                    _ => 24.0,
                };
                return RiseTransitSet {
                    has_transit: true,
                    has_rise: has_r,
                    has_set: has_s,
                    jd_rise_ut: opt_r,
                    jd_transit_ut: Some(t_trnsit),
                    jd_set_ut: opt_s,
                    rise_hms_ut: r_hms,
                    transit_hms_ut: jd_to_hms_ut(t_trnsit),
                    set_hms_ut: s_hms,
                    rise_cal_ut: r_cal,
                    transit_cal_ut: tr_cal,
                    set_cal_ut: s_cal,
                    transit_elevation_deg: elev_trnsit,
                    visible_hours: vis,
                };
            }
        }
        return RiseTransitSet {
            has_transit: true,
            has_rise: false,
            has_set: false,
            jd_rise_ut: None,
            jd_transit_ut: Some(t_trnsit),
            jd_set_ut: None,
            rise_hms_ut: "Circumpolar / Não nasce".to_string(),
            transit_hms_ut: jd_to_hms_ut(t_trnsit),
            set_hms_ut: "Circumpolar / Não se põe".to_string(),
            rise_cal_ut: "Circumpolar".to_string(),
            transit_cal_ut: tr_cal,
            set_cal_ut: "Circumpolar".to_string(),
            transit_elevation_deg: elev_trnsit,
            visible_hours: if elev_trnsit > 0.0 { 24.0 } else { 0.0 },
        };
    }

    // 2. Itera o horário de Nascer (Rise) anterior ao trânsito (`iter_trnsit` linhas 224-262)
    let mut date_r = date_trnsit;
    let mut t_r = date_r + r_rise_init / tpi;
    if t_r >= t_trnsit {
        date_r -= 1.0;
        t_r = date_r + r_rise_init / tpi;
    }
    let mut has_rise = true;
    for _ in 0..11 {
        let t0 = t_r;
        let (ok, _, r_ri, _, _, _, _) = eval_at_ut(t0);
        if !ok {
            has_rise = false;
            break;
        }
        t_r = date_r + r_ri / tpi;
        if t_r > t_trnsit {
            date_r -= 1.0;
            t_r = date_r + r_ri / tpi;
        }
        if (t_r - t0).abs() <= 0.0001 {
            break;
        }
    }

    // 3. Itera o horário de Pôr (Set) posterior ao trânsito (`iter_trnsit` linhas 265-302)
    let mut date_s = date_trnsit;
    let mut t_s = date_s + r_set_init / tpi;
    if t_s <= t_trnsit {
        date_s += 1.0;
        t_s = date_s + r_set_init / tpi;
    }
    let mut has_set = true;
    for _ in 0..11 {
        let t0 = t_s;
        let (ok, _, _, r_se, _, _, _) = eval_at_ut(t0);
        if !ok {
            has_set = false;
            break;
        }
        t_s = date_s + r_se / tpi;
        if t_s < t_trnsit {
            date_s += 1.0;
            t_s = date_s + r_se / tpi;
        }
        if (t_s - t0).abs() <= 0.0001 {
            break;
        }
    }

    if !has_rise || !has_set {
        let mut eval_elev = |t: f64| -> f64 { eval_at_ut(t).5 };
        let (opt_r, opt_s) = no_rise_set(t_trnsit, elev_true_trnsit, elev_thresh, &mut eval_elev);
        if let Some(r_val) = opt_r {
            has_rise = true;
            t_r = r_val;
        }
        if let Some(s_val) = opt_s {
            has_set = true;
            t_s = s_val;
        }
    }

    let visible_hours = if has_rise && has_set && t_s > t_r && (t_s - t_r) < 1.0 {
        24.0 * (t_s - t_r)
    } else if elev_trnsit > 0.0 {
        24.0
    } else {
        0.0
    };

    let rise_cal = if has_rise {
        jtocal_parts(t_r).6
    } else {
        "—".to_string()
    };
    let set_cal = if has_set {
        jtocal_parts(t_s).6
    } else {
        "—".to_string()
    };

    RiseTransitSet {
        has_transit: true,
        has_rise,
        has_set,
        jd_rise_ut: if has_rise { Some(t_r) } else { None },
        jd_transit_ut: Some(t_trnsit),
        jd_set_ut: if has_set { Some(t_s) } else { None },
        rise_hms_ut: if has_rise {
            jd_to_hms_ut(t_r)
        } else {
            "—".to_string()
        },
        transit_hms_ut: jd_to_hms_ut(t_trnsit),
        set_hms_ut: if has_set {
            jd_to_hms_ut(t_s)
        } else {
            "—".to_string()
        },
        rise_cal_ut: rise_cal,
        transit_cal_ut: tr_cal,
        set_cal_ut: set_cal,
        transit_elevation_deg: elev_trnsit,
        visible_hours,
    }
}

