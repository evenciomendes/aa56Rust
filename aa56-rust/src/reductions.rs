//! Astrometric and Apparent Reductions (`sun.c`, `rplanet.c`, `lightt.c`, `deflec.c`, `annuab.c`, `angles.c`).
//! Combines all DE404 planetary ephemerides, `domoon` lunar ephemeris, Lunar Nodes / Black Moon Lilith,
//! `2060 Chiron`, `conjunct` lunar phases/seasons, and `moonrise` preview into `EphemerisReport`.

use serde::Serialize;

use crate::conjunct::{compute_conjunct_report, ConjunctReport};
use crate::constel::whatconstel;
use crate::domoon::{domoon, DomoonResult};
use crate::earth_orientation::{
    epsiln, geodetic_to_geocentric, lonlat, nutate, nutlo, precess, sidrlt, topocentric_altaz,
};
use crate::kepler::{
    builtin_planets, compute_lunar_nodes_and_lilith, kepler, velearth, De404Body, OrbitElement,
};
use crate::math_time::{
    caltoj, deltat, format_dms, format_hms_rad, jtocal_parts, jtocal_string, modtp, zatan2,
    AEARTH_M, AU_KM, CLIGHT_AUD, DTR, J2000, RTD, RTS, STR, TPI,
};
use crate::moonrise::{compute_moonrise_table, MoonriseReport};
use crate::transit::{
    compute_rise_transit_set, RiseSetTransit, HORIZON_PLANET_STAR, HORIZON_SUN,
};

#[derive(Debug, Clone)]
pub struct EphemerisInput {
    pub year: i32,
    pub month: i32,
    pub day: f64,
    pub hour: f64,
    pub minute: f64,
    pub second: f64,
    pub delta_t_override: Option<f64>,
    pub obs_lon_deg: f64,
    pub obs_lat_deg: f64,
    pub obs_height_m: f64,
    pub temp_c: f64,
    pub pressure_mb: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct AnglesInfo {
    pub r_helio_au: f64,
    pub r_geo_au: f64,
    pub r_earth_au: f64,
    pub phase_factor: f64,
    pub elongation_deg: f64,
    pub phase_angle_deg: f64,
    pub light_time_min: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct BodyEphemeris {
    pub id: i32,
    pub name: String,
    pub category: String,
    pub ephem_source: String,
    pub constellation: String,
    pub helio_rect_j2000: [f64; 3],
    pub helio_lon_deg: f64,
    pub helio_lat_deg: f64,
    pub helio_dist_au: f64,
    pub geo_dist_au: f64,
    pub light_time_min: f64,
    pub astrometric_j2000_ra_rad: f64,
    pub astrometric_j2000_dec_rad: f64,
    pub astrometric_j2000_ra_hms: String,
    pub astrometric_j2000_dec_dms: String,
    pub apparent_ra_rad: f64,
    pub apparent_dec_rad: f64,
    pub apparent_ra_hms: String,
    pub apparent_dec_dms: String,
    pub apparent_ecl_lon_deg: f64,
    pub apparent_ecl_lat_deg: f64,
    pub apparent_ecl_lon_dms: String,
    pub apparent_ecl_lat_dms: String,
    pub zodiac_sign: String,
    pub zodiac_deg_dms: String,
    pub horiz_parallax_arcsec: f64,
    pub semidiameter_arcsec: f64,
    pub visual_magnitude: f64,
    pub illuminated_fraction: f64,
    pub solar_elongation_deg: f64,
    pub phase_note: String,
    pub topo_ra_hms: String,
    pub topo_dec_dms: String,
    pub topo_alt_deg: f64,
    pub topo_az_deg: f64,
    pub refracted_alt_deg: f64,
    pub lha_hours: f64,
    pub rise_transit_set: RiseSetTransit,
}

#[derive(Debug, Clone, Serialize)]
pub struct EphemerisReport {
    pub engine: String,
    pub jd_ut: f64,
    pub jd_tdt: f64,
    pub delta_t_sec: f64,
    pub ut_date_str: String,
    pub tdt_date_str: String,
    pub gst_hours: f64,
    pub gst_hms: String,
    pub lst_hours: f64,
    pub lst_hms: String,
    pub mean_obliquity_deg: f64,
    pub mean_obliquity_dms: String,
    pub nutation_lon_arcsec: f64,
    pub nutation_obliq_arcsec: f64,
    pub obs_lon_deg: f64,
    pub obs_lat_deg: f64,
    pub obs_height_m: f64,
    pub temp_c: f64,
    pub pressure_mb: f64,
    pub earth_helio_j2000: [f64; 3],
    pub earth_helio_lon_deg: f64,
    pub earth_helio_lat_deg: f64,
    pub earth_helio_dist_au: f64,
    pub bodies: Vec<BodyEphemeris>,
    pub domoon_details: DomoonResult,
    pub conjunct_summary: ConjunctReport,
    pub moonrise_preview: MoonriseReport,
}

/// Zodiac sign helper (`0..12` from ecliptic longitude in degrees).
pub fn format_zodiac(lon_deg: f64) -> (String, String) {
    let signs = [
        "Áries ♈",
        "Touro ♉",
        "Gêmeos ♊",
        "Câncer ♋",
        "Leão ♌",
        "Virgem ♍",
        "Libra ♎",
        "Escorpião ♏",
        "Sagitário ♐",
        "Capricórnio ♑",
        "Aquário ♒",
        "Peixes ♓",
    ];
    let norm = ((lon_deg % 360.0) + 360.0) % 360.0;
    let idx = ((norm / 30.0).floor() as usize) % 12;
    let rem_deg = norm - (idx as f64) * 30.0;
    let deg = rem_deg.floor() as i32;
    let min_f = (rem_deg - (deg as f64)) * 60.0;
    let min = min_f.floor() as i32;
    let sec = (min_f - (min as f64)) * 60.0;
    (
        signs[idx].to_string(),
        format!("{:02}°{:02}'{:04.1}\"", deg, min, sec),
    )
}

/// Annual aberration (`annuab.c`) using `velearth(jed)`.
/// Modifies unit direction vector `p` (J2000 equatorial rectangular) in place and returns displacement `dp`.
pub fn annuab(p: &mut [f64; 3], jed: f64) -> [f64; 3] {
    let vearth = velearth(jed);
    let mut v = [0.0f64; 3];
    let mut betai = 0.0f64;
    let mut p_dot_v = 0.0f64;
    for i in 0..3 {
        let vi = vearth[i] / CLIGHT_AUD;
        v[i] = vi;
        betai += vi * vi;
        p_dot_v += p[i] * vi;
    }
    betai = (1.0 - betai).sqrt();
    let c = 1.0 + p_dot_v;
    let a = betai / c;
    let b = (1.0 + p_dot_v / (1.0 + betai)) / c;

    let mut r = 0.0f64;
    let mut dp = [0.0f64; 3];
    for i in 0..3 {
        let x = a * p[i] + b * v[i];
        dp[i] = x - p[i];
        p[i] = x;
        r += x * x;
    }
    r = r.sqrt();
    for i in 0..3 {
        p[i] /= r;
    }
    dp
}

/// Relativistic deflection of light due to the Sun's gravitational field (`deflec.c`).
/// `p`: geocentric direction vector (unit), `q`: heliocentric object vector, `e`: heliocentric Earth vector.
pub fn relativity(p: &mut [f64; 3], q: &[f64; 3], e: &[f64; 3]) -> [f64; 3] {
    let mut pq = 0.0;
    let mut ep = 0.0;
    let mut qe = 0.0;

    let mut e_norm = 0.0;
    let mut q_norm = 0.0;
    let mut p_norm = 0.0;
    for i in 0..3 {
        e_norm += e[i] * e[i];
        q_norm += q[i] * q[i];
        p_norm += p[i] * p[i];
    }
    let eo = e_norm.sqrt();
    let qo = q_norm.sqrt();
    let po = p_norm.sqrt();

    if eo == 0.0 || qo == 0.0 || po == 0.0 {
        return [0.0; 3];
    }

    for i in 0..3 {
        ep += (e[i] / eo) * (p[i] / po);
        pq += (p[i] / po) * (q[i] / qo);
        qe += (q[i] / qo) * (e[i] / eo);
    }

    if (1.0 + qe).abs() < 1.0e-5 {
        return [0.0; 3];
    }

    let c = 1.974e-8 / (eo * (1.0 + qe));
    let mut dp = [0.0f64; 3];
    let mut r = 0.0;
    for i in 0..3 {
        dp[i] = c * (pq * (e[i] / eo) - ep * (q[i] / qo));
        p[i] = (p[i] / po) + dp[i];
        r += p[i] * p[i];
    }
    r = r.sqrt();
    for i in 0..3 {
        p[i] /= r;
    }
    dp
}

/// Sun-Earth-Object Euclidean distances and cosine angles (`angles.c`).
pub struct AnglesData {
    pub r_helio_au: f64,
    pub r_geo_au: f64,
    pub r_earth_au: f64,
    pub ep_cos: f64,
    pub pq_cos: f64,
    pub qe_cos: f64,
    pub phase_factor: f64,
}

pub fn angles(p: &[f64; 3], q: &[f64; 3], e: &[f64; 3]) -> AnglesData {
    let mut a = 0.0;
    let mut b = 0.0;
    let mut c = 0.0;
    let mut ep = 0.0;
    let mut pq = 0.0;
    let mut qe = 0.0;
    for i in 0..3 {
        a += e[i] * e[i];
        b += p[i] * p[i];
        c += q[i] * q[i];
        ep += e[i] * p[i];
        pq += p[i] * q[i];
        qe += q[i] * e[i];
    }
    let eo = a.sqrt();
    let po = b.sqrt();
    let qo = c.sqrt();

    let ep_cos = if eo * po > 0.0 { ep / (eo * po) } else { 0.0 };
    let pq_cos = if po * qo > 0.0 { pq / (po * qo) } else { 1.0 };
    let qe_cos = if qo * eo > 0.0 { qe / (qo * eo) } else { 0.0 };
    let phase_factor = 0.5 * (1.0 + pq_cos);

    AnglesData {
        r_helio_au: qo,
        r_geo_au: po,
        r_earth_au: eo,
        ep_cos,
        pq_cos,
        qe_cos,
        phase_factor,
    }
}

/// Light-time retardation for planetary / orbital body (`lightt.c`).
/// Returns `(p_unit_j2000, q_retarded_j2000, p_vec_j2000, angles_data, light_time_days)`.
pub fn lightt(
    elem: &mut OrbitElement,
    earth_j2000: &[f64; 3],
    jed: f64,
) -> ([f64; 3], [f64; 3], [f64; 3], AnglesData, f64) {
    let mut p = [0.0f64; 3];
    let mut p0 = [0.0f64; 3];
    let mut ptemp = [0.0f64; 3];
    let mut q = [0.0f64; 3];
    let t_light;

    // Save initial geometric position at `jed`
    let (q0, _) = kepler(jed, elem);
    for i in 0..3 {
        p0[i] = q0[i] - earth_j2000[i];
    }

    let vearth = velearth(jed);

    for _ in 0..2 {
        let mut a = 0.0;
        for i in 0..3 {
            let x = q0[i] - earth_j2000[i];
            p[i] = x;
            a += x * x;
        }
        t_light = a.sqrt() / CLIGHT_AUD;
        let (q_ret, _) = kepler(jed - t_light, elem);
        q = q_ret;
    }

    // Final light-time and aberration correction in `lightt.c`
    let mut a = 0.0;
    for i in 0..3 {
        let x = q[i] - earth_j2000[i];
        p[i] = x;
        a += x * x;
    }
    t_light = a.sqrt() / CLIGHT_AUD;

    for i in 0..3 {
        ptemp[i] = q[i] - (earth_j2000[i] - t_light * vearth[i]);
        p[i] = q[i] - earth_j2000[i];
    }

    let ang = angles(&p, &q0, earth_j2000);

    let mut r = 0.0;
    for i in 0..3 {
        r += ptemp[i] * ptemp[i];
    }
    r = r.sqrt();
    let mut p_unit = [0.0f64; 3];
    if r > 0.0 {
        for i in 0..3 {
            p_unit[i] = ptemp[i] / r;
        }
    }

    (p_unit, q, p, ang, t_light)
}

/// Reduce a planetary / Keplerian body at `jd_tdt` (`rplanet.c` -> `reduce` + `doplanet`).
pub fn reduce_planet(
    id: i32,
    category: &str,
    elem: &mut OrbitElement,
    earth_j2000: &[f64; 3],
    jd_tdt: f64,
    jd_ut: f64,
    obs_lon_deg: f64,
    obs_lat_deg: f64,
    obs_height_m: f64,
    temp_c: f64,
    pressure_mb: f64,
) -> BodyEphemeris {
    let ephem_source = if elem.de404 != De404Body::None {
        "DE404 (gplan.c)".to_string()
    } else {
        "Keplerian (kepler.c)".to_string()
    };

    // Geometric heliocentric coordinates at TDT
    let (helio_rect, helio_pol) = kepler(jd_tdt, elem);

    // Light time and aberration (`lightt.c`)
    let (mut p, q, _p_vec, ang, t_light) = lightt(elem, earth_j2000, jd_tdt);

    // Astrometric J2000 RA/Dec (after light time, before gravitational deflection/aberration)
    let mut p_astro = [
        q[0] - earth_j2000[0],
        q[1] - earth_j2000[1],
        q[2] - earth_j2000[2],
    ];
    let astro_norm =
        (p_astro[0] * p_astro[0] + p_astro[1] * p_astro[1] + p_astro[2] * p_astro[2]).sqrt();
    if astro_norm > 0.0 {
        for i in 0..3 {
            p_astro[i] /= astro_norm;
        }
    }
    let astro_ra = zatan2(p_astro[0], p_astro[1]);
    let astro_dec = p_astro[2].clamp(-1.0, 1.0).asin();
    let constel = whatconstel(&p_astro, J2000).to_string();

    // Relativistic light deflection (`deflec.c`)
    relativity(&mut p, &q, earth_j2000);

    // Note: `lightt` already applied `-t * vearth[i]` (annual + planetary aberration together),
    // which is Moshier's exact reduction in `rplanet.c`!
    // Precess from J2000 to equinox and ecliptic of date (`precess(p, TDT, -1)`)
    precess(&mut p, jd_tdt, -1);

    // Ecliptic longitude/latitude of date (`lonlat`) before nutation, plus nutation in lon
    let (ecl_lon_mean, ecl_lat_mean, _) = lonlat(&p, jd_tdt);
    let (nutl, _nuto) = nutlo(jd_tdt);
    let app_ecl_lon_rad = modtp(ecl_lon_mean + nutl);
    let app_ecl_lat_rad = ecl_lat_mean;

    // Nutation (`nutate(TDT, p)`)
    nutate(jd_tdt, &mut p);

    let app_ra = zatan2(p[0], p[1]);
    let app_dec = p[2].clamp(-1.0, 1.0).asin();

    // Horizontal parallax (`8.794 / (RTS * robject)` in `rplanet.c`) and Semidiameter
    let hp_rad = if ang.r_geo_au > 0.0 {
        8.794 / (RTS * ang.r_geo_au)
    } else {
        0.0
    };
    let sd_rad = if ang.r_geo_au > 0.0 {
        elem.sdiam / (RTS * ang.r_geo_au)
    } else {
        0.0
    };

    // Visual magnitude (`rplanet.c` line 139: `elem->mag + 5.0*log10(r*robject) - 2.5*log10(fac)`)
    let a_prod = ang.r_helio_au * ang.r_geo_au;
    let mut mag = elem.mag;
    if a_prod > 0.0 {
        mag += 5.0 * a_prod.log10();
    }
    if ang.phase_factor > 1.0e-6 {
        mag -= 2.5 * ang.phase_factor.log10();
    }

    let solar_elong_deg = RTD * (-ang.ep_cos).clamp(-1.0, 1.0).acos();
    let (zodiac_sign, zodiac_deg_dms) = format_zodiac(app_ecl_lon_rad * RTD);

    let geo_obs = geodetic_to_geocentric(obs_lat_deg, obs_height_m);
    let topo = topocentric_altaz(
        app_ra,
        app_dec,
        jd_ut,
        jd_tdt,
        obs_lon_deg,
        obs_lat_deg,
        &geo_obs,
        hp_rad,
        temp_c,
        pressure_mb,
    );

    // Rise / Transit / Set (`trnsit.c`)
    let mut elem_clone = elem.clone();
    let mut earth_el = OrbitElement {
        obname: "Earth".to_string(),
        epoch: 2446800.5,
        i: 0.0,
        W: 0.0,
        w: 102.884,
        a: 0.999999,
        dm: 0.985611,
        ecc: 0.016713,
        M: 1.1791,
        equinox: 2446800.5,
        mag: -4.0,
        sdiam: 8.79,
        L: 0.0,
        r: 0.0,
        plat: 0.0,
        de404: De404Body::EarthBary,
    };
    let dt_days = jd_tdt - jd_ut;
    let rts = compute_rise_transit_set(
        jd_ut,
        obs_lon_deg,
        obs_lat_deg,
        &geo_obs,
        HORIZON_PLANET_STAR,
        temp_c,
        pressure_mb,
        |t_ut| {
            let t_tdt = t_ut + dt_days;
            let (re, _) = kepler(t_tdt, &mut earth_el);
            let (mut pu, qu, _, ang_i, _) = lightt(&mut elem_clone, &re, t_tdt);
            relativity(&mut pu, &qu, &re);
            precess(&mut pu, t_tdt, -1);
            nutate(t_tdt, &mut pu);
            let r_ra = zatan2(pu[0], pu[1]);
            let r_dec = pu[2].clamp(-1.0, 1.0).asin();
            let r_hp = if ang_i.r_geo_au > 0.0 {
                8.794 / (RTS * ang_i.r_geo_au)
            } else {
                0.0
            };
            let r_sd = if ang_i.r_geo_au > 0.0 {
                elem_clone.sdiam / (RTS * ang_i.r_geo_au)
            } else {
                0.0
            };
            (r_ra, r_dec, r_hp, r_sd)
        },
    );

    BodyEphemeris {
        id,
        name: elem.obname.clone(),
        category: category.to_string(),
        ephem_source,
        constellation: constel,
        helio_rect_j2000: helio_rect,
        helio_lon_deg: helio_pol[0] * RTD,
        helio_lat_deg: helio_pol[1] * RTD,
        helio_dist_au: helio_pol[2],
        geo_dist_au: ang.r_geo_au,
        light_time_min: t_light * 1440.0,
        astrometric_j2000_ra_rad: astro_ra,
        astrometric_j2000_dec_rad: astro_dec,
        astrometric_j2000_ra_hms: format_hms_rad(astro_ra),
        astrometric_j2000_dec_dms: format_dms(astro_dec),
        apparent_ra_rad: app_ra,
        apparent_dec_rad: app_dec,
        apparent_ra_hms: format_hms_rad(app_ra),
        apparent_dec_dms: format_dms(app_dec),
        apparent_ecl_lon_deg: app_ecl_lon_rad * RTD,
        apparent_ecl_lat_deg: app_ecl_lat_rad * RTD,
        apparent_ecl_lon_dms: format_dms(app_ecl_lon_rad),
        apparent_ecl_lat_dms: format_dms(app_ecl_lat_rad),
        zodiac_sign,
        zodiac_deg_dms,
        horiz_parallax_arcsec: hp_rad * RTS,
        semidiameter_arcsec: sd_rad * RTS,
        visual_magnitude: mag,
        illuminated_fraction: ang.phase_factor,
        solar_elongation_deg: solar_elong_deg,
        phase_note: format!("Ilum: {:.1}%", ang.phase_factor * 100.0),
        topo_ra_hms: format_hms_rad(topo.topo_ra_rad),
        topo_dec_dms: format_dms(topo.topo_dec_rad),
        topo_alt_deg: topo.true_alt_deg,
        topo_az_deg: topo.true_az_deg,
        refracted_alt_deg: topo.refracted_alt_deg,
        lha_hours: topo.lha_hours,
        rise_transit_set: rts,
    }
}

/// Compute Sun apparent ephemeris (`dosun` in `sun.c`).
pub fn reduce_sun(
    earth_j2000: &[f64; 3],
    earth_polar_j2000: &[f64; 3],
    jd_tdt: f64,
    jd_ut: f64,
    obs_lon_deg: f64,
    obs_lat_deg: f64,
    obs_height_m: f64,
    temp_c: f64,
    pressure_mb: f64,
) -> BodyEphemeris {
    let mut earth_el = OrbitElement {
        obname: "Earth".to_string(),
        epoch: 2446800.5,
        i: 0.0,
        W: 0.0,
        w: 102.884,
        a: 0.999999,
        dm: 0.985611,
        ecc: 0.016713,
        M: 1.1791,
        equinox: 2446800.5,
        mag: -4.0,
        sdiam: 8.79,
        L: 0.0,
        r: 0.0,
        plat: 0.0,
        de404: De404Body::EarthBary,
    };

    // Light-time iteration in `sun.c` (`t = 0.0057755183 * r`)
    let mut r = earth_polar_j2000[2];
    let mut t_light = 0.0;
    let mut ecr = *earth_j2000;
    for _ in 0..2 {
        t_light = 0.0057755183 * r;
        let (ecr_ret, pol_ret) = kepler(jd_tdt - t_light, &mut earth_el);
        ecr = ecr_ret;
        r = pol_ret[2];
    }

    let mut s_vec = [-ecr[0], -ecr[1], -ecr[2]];
    let s_norm = (s_vec[0] * s_vec[0] + s_vec[1] * s_vec[1] + s_vec[2] * s_vec[2]).sqrt();
    let mut p_astro = [s_vec[0] / s_norm, s_vec[1] / s_norm, s_vec[2] / s_norm];
    let astro_ra = zatan2(p_astro[0], p_astro[1]);
    let astro_dec = p_astro[2].clamp(-1.0, 1.0).asin();
    let constel = whatconstel(&p_astro, J2000).to_string();

    // Annual aberration (`annuab`)
    annuab(&mut p_astro, jd_tdt);
    // Precess to date (`precess(s, TDT, -1)`)
    precess(&mut p_astro, jd_tdt, -1);

    let (ecl_lon_mean, ecl_lat_mean, _) = lonlat(&p_astro, jd_tdt);
    let (nutl, _nuto) = nutlo(jd_tdt);
    let app_ecl_lon_rad = modtp(ecl_lon_mean + nutl);
    let app_ecl_lat_rad = ecl_lat_mean;

    // Nutation (`nutate(TDT, s)`)
    nutate(jd_tdt, &mut p_astro);

    let app_ra = zatan2(p_astro[0], p_astro[1]);
    let app_dec = p_astro[2].clamp(-1.0, 1.0).asin();

    let hp_rad = 8.794 / (RTS * r);
    let sd_rad = 959.63 / (RTS * r);
    let (zodiac_sign, zodiac_deg_dms) = format_zodiac(app_ecl_lon_rad * RTD);

    let geo_obs = geodetic_to_geocentric(obs_lat_deg, obs_height_m);
    let topo = topocentric_altaz(
        app_ra,
        app_dec,
        jd_ut,
        jd_tdt,
        obs_lon_deg,
        obs_lat_deg,
        &geo_obs,
        hp_rad,
        temp_c,
        pressure_mb,
    );

    let dt_days = jd_tdt - jd_ut;
    let rts = compute_rise_transit_set(
        jd_ut,
        obs_lon_deg,
        obs_lat_deg,
        &geo_obs,
        HORIZON_SUN,
        temp_c,
        pressure_mb,
        |t_ut| {
            let t_tdt = t_ut + dt_days;
            let (_e0, pol0) = kepler(t_tdt, &mut earth_el);
            let tl = 0.0057755183 * pol0[2];
            let (ec_r, pol_r) = kepler(t_tdt - tl, &mut earth_el);
            let mut sv = [-ec_r[0], -ec_r[1], -ec_r[2]];
            let sn = (sv[0] * sv[0] + sv[1] * sv[1] + sv[2] * sv[2]).sqrt();
            for k in 0..3 {
                sv[k] /= sn;
            }
            annuab(&mut sv, t_tdt);
            precess(&mut sv, t_tdt, -1);
            nutate(t_tdt, &mut sv);
            let r_ra = zatan2(sv[0], sv[1]);
            let r_dec = sv[2].clamp(-1.0, 1.0).asin();
            let r_hp = 8.794 / (RTS * pol_r[2]);
            let r_sd = 959.63 / (RTS * pol_r[2]);
            (r_ra, r_dec, r_hp, r_sd)
        },
    );

    BodyEphemeris {
        id: 0,
        name: "Sun".to_string(),
        category: "Luminary".to_string(),
        ephem_source: "DE404 (ear404 + embofs + sun.c)".to_string(),
        constellation: constel,
        helio_rect_j2000: [0.0, 0.0, 0.0],
        helio_lon_deg: earth_polar_j2000[0] * RTD,
        helio_lat_deg: earth_polar_j2000[1] * RTD,
        helio_dist_au: earth_polar_j2000[2],
        geo_dist_au: r,
        light_time_min: t_light * 1440.0,
        astrometric_j2000_ra_rad: astro_ra,
        astrometric_j2000_dec_rad: astro_dec,
        astrometric_j2000_ra_hms: format_hms_rad(astro_ra),
        astrometric_j2000_dec_dms: format_dms(astro_dec),
        apparent_ra_rad: app_ra,
        apparent_dec_rad: app_dec,
        apparent_ra_hms: format_hms_rad(app_ra),
        apparent_dec_dms: format_dms(app_dec),
        apparent_ecl_lon_deg: app_ecl_lon_rad * RTD,
        apparent_ecl_lat_deg: app_ecl_lat_rad * RTD,
        apparent_ecl_lon_dms: format_dms(app_ecl_lon_rad),
        apparent_ecl_lat_dms: format_dms(app_ecl_lat_rad),
        zodiac_sign,
        zodiac_deg_dms,
        horiz_parallax_arcsec: hp_rad * RTS,
        semidiameter_arcsec: sd_rad * RTS,
        visual_magnitude: -26.74,
        illuminated_fraction: 1.0,
        solar_elongation_deg: 0.0,
        phase_note: "Solar Disk".to_string(),
        topo_ra_hms: format_hms_rad(topo.topo_ra_rad),
        topo_dec_dms: format_dms(topo.topo_dec_rad),
        topo_alt_deg: topo.true_alt_deg,
        topo_az_deg: topo.true_az_deg,
        refracted_alt_deg: topo.refracted_alt_deg,
        lha_hours: topo.lha_hours,
        rise_transit_set: rts,
    }
}

/// Convert ecliptic longitude of date (with `lat = 0`) to a mathematical point `BodyEphemeris`
/// (used for Lunar Nodes & Black Moon Lilith).
fn make_ecliptic_point_body(
    id: i32,
    name: &str,
    source: &str,
    lon_rad_of_date: f64,
    jd_tdt: f64,
    jd_ut: f64,
    obs_lon_deg: f64,
    obs_lat_deg: f64,
    obs_height_m: f64,
    temp_c: f64,
    pressure_mb: f64,
) -> BodyEphemeris {
    let eps = epsiln(jd_tdt);
    let mut p = [
        lon_rad_of_date.cos(),
        lon_rad_of_date.sin() * eps.cos(),
        lon_rad_of_date.sin() * eps.sin(),
    ];

    let mut p_j2000 = p;
    precess(&mut p_j2000, jd_tdt, 1);
    let astro_ra = zatan2(p_j2000[0], p_j2000[1]);
    let astro_dec = p_j2000[2].clamp(-1.0, 1.0).asin();
    let constel = whatconstel(&p_j2000, J2000).to_string();

    let (nutl, _) = nutlo(jd_tdt);
    let app_lon = modtp(lon_rad_of_date + nutl);
    nutate(jd_tdt, &mut p);
    let app_ra = zatan2(p[0], p[1]);
    let app_dec = p[2].clamp(-1.0, 1.0).asin();
    let (zodiac_sign, zodiac_deg_dms) = format_zodiac(app_lon * RTD);

    let geo_obs = geodetic_to_geocentric(obs_lat_deg, obs_height_m);
    let topo = topocentric_altaz(
        app_ra,
        app_dec,
        jd_ut,
        jd_tdt,
        obs_lon_deg,
        obs_lat_deg,
        &geo_obs,
        0.0,
        temp_c,
        pressure_mb,
    );

    let rts = compute_rise_transit_set(
        jd_ut,
        obs_lon_deg,
        obs_lat_deg,
        &geo_obs,
        HORIZON_PLANET_STAR,
        temp_c,
        pressure_mb,
        |_t_ut| (app_ra, app_dec, 0.0, 0.0),
    );

    BodyEphemeris {
        id,
        name: name.to_string(),
        category: "Lunar Point".to_string(),
        ephem_source: source.to_string(),
        constellation: constel,
        helio_rect_j2000: [0.0, 0.0, 0.0],
        helio_lon_deg: app_lon * RTD,
        helio_lat_deg: 0.0,
        helio_dist_au: 0.0025695,
        geo_dist_au: 0.0025695,
        light_time_min: 0.0,
        astrometric_j2000_ra_rad: astro_ra,
        astrometric_j2000_dec_rad: astro_dec,
        astrometric_j2000_ra_hms: format_hms_rad(astro_ra),
        astrometric_j2000_dec_dms: format_dms(astro_dec),
        apparent_ra_rad: app_ra,
        apparent_dec_rad: app_dec,
        apparent_ra_hms: format_hms_rad(app_ra),
        apparent_dec_dms: format_dms(app_dec),
        apparent_ecl_lon_deg: app_lon * RTD,
        apparent_ecl_lat_deg: 0.0,
        apparent_ecl_lon_dms: format_dms(app_lon),
        apparent_ecl_lat_dms: format_dms(0.0),
        zodiac_sign,
        zodiac_deg_dms,
        horiz_parallax_arcsec: 0.0,
        semidiameter_arcsec: 0.0,
        visual_magnitude: 0.0,
        illuminated_fraction: 0.0,
        solar_elongation_deg: 0.0,
        phase_note: "Ecliptic Point".to_string(),
        topo_ra_hms: format_hms_rad(topo.topo_ra_rad),
        topo_dec_dms: format_dms(topo.topo_dec_rad),
        topo_alt_deg: topo.true_alt_deg,
        topo_az_deg: topo.true_az_deg,
        refracted_alt_deg: topo.refracted_alt_deg,
        lha_hours: topo.lha_hours,
        rise_transit_set: rts,
    }
}

/// Complete Ephemeris Report for `aa-56` (`aa.c` + `domoon.c` + `conjunct.c` + `moonrise.c`).
pub fn compute_ephemeris_report(input: &EphemerisInput) -> EphemerisReport {
    let day_frac =
        input.day + (input.hour + input.minute / 60.0 + input.second / 3600.0) / 24.0;
    let jd_ut = caltoj(input.year, input.month, day_frac);

    let y_dec = (input.year as f64)
        + ((input.month as f64) - 1.0) / 12.0
        + day_frac / 365.25;
    let delta_t_sec = input.delta_t_override.unwrap_or_else(|| deltat(y_dec));
    let jd_tdt = jd_ut + delta_t_sec / 86400.0;

    // Earth heliocentric J2000 position (`ear404` + `embofs`)
    let mut planets = builtin_planets();
    let mut earth_el = planets[2].clone();
    let (earth_helio_j2000, earth_polar_j2000) = kepler(jd_tdt, &mut earth_el);

    let gst_rad = sidrlt(jd_ut, jd_tdt);
    let lst_rad = modtp(gst_rad + input.obs_lon_deg * DTR);
    let eps = epsiln(jd_tdt);
    let (nutl, nuto) = nutlo(jd_tdt);

    let mut bodies = Vec::new();

    // 0. Sun (`sun.c`)
    bodies.push(reduce_sun(
        &earth_helio_j2000,
        &earth_polar_j2000,
        jd_tdt,
        jd_ut,
        input.obs_lon_deg,
        input.obs_lat_deg,
        input.obs_height_m,
        input.temp_c,
        input.pressure_mb,
    ));

    // 1..9 Planets (`Mercury`, `Venus`, `Mars`, `Jupiter`, `Saturn`, `Uranus`, `Neptune`, `Pluto`)
    let mut body_id = 1;
    for p in planets.iter_mut() {
        if p.de404 == De404Body::EarthBary {
            continue;
        }
        let cat = if p.de404 != De404Body::None {
            "Major Planet"
        } else {
            "Centaur / Minor Body"
        };
        bodies.push(reduce_planet(
            body_id,
            cat,
            p,
            &earth_helio_j2000,
            jd_tdt,
            jd_ut,
            input.obs_lon_deg,
            input.obs_lat_deg,
            input.obs_height_m,
            input.temp_c,
            input.pressure_mb,
        ));
        body_id += 1;
    }

    // Moon (`domoon.c` with `mlr404` + `mlat404`)
    let dm = domoon(
        jd_tdt,
        jd_ut,
        input.obs_lon_deg,
        input.obs_lat_deg,
        input.obs_height_m,
        input.temp_c,
        input.pressure_mb,
    );
    let mut moon_j2000 = dm.geo_rect_of_date;
    precess(&mut moon_j2000, jd_tdt, 1);
    let moon_constel = whatconstel(&moon_j2000, J2000).to_string();
    let (moon_zodiac, moon_zodiac_dms) = format_zodiac(dm.app_ecl_lon_deg);
    let moon_mag = -12.74
        + 0.026 * (180.0 - dm.solar_elongation_deg).abs()
        + 4.0e-9 * (180.0 - dm.solar_elongation_deg).powi(4);

    bodies.insert(
        1,
        BodyEphemeris {
            id: 10,
            name: "Moon".to_string(),
            category: "Luminary".to_string(),
            ephem_source: "DE404 Lunar (mlr404 + mlat404 + domoon.c)".to_string(),
            constellation: moon_constel,
            helio_rect_j2000: moon_j2000,
            helio_lon_deg: dm.app_ecl_lon_deg,
            helio_lat_deg: dm.app_ecl_lat_deg,
            helio_dist_au: dm.distance_au,
            geo_dist_au: dm.distance_au,
            light_time_min: (dm.distance_au / CLIGHT_AUD) * 1440.0,
            astrometric_j2000_ra_rad: zatan2(moon_j2000[0], moon_j2000[1]),
            astrometric_j2000_dec_rad: moon_j2000[2].clamp(-1.0, 1.0).asin(),
            astrometric_j2000_ra_hms: format_hms_rad(zatan2(moon_j2000[0], moon_j2000[1])),
            astrometric_j2000_dec_dms: format_dms(moon_j2000[2].clamp(-1.0, 1.0).asin()),
            apparent_ra_rad: dm.app_ra_rad,
            apparent_dec_rad: dm.app_dec_rad,
            apparent_ra_hms: dm.app_ra_hms.clone(),
            apparent_dec_dms: dm.app_dec_dms.clone(),
            apparent_ecl_lon_deg: dm.app_ecl_lon_deg,
            apparent_ecl_lat_deg: dm.app_ecl_lat_deg,
            apparent_ecl_lon_dms: dm.app_ecl_lon_dms.clone(),
            apparent_ecl_lat_dms: dm.app_ecl_lat_dms.clone(),
            zodiac_sign: moon_zodiac,
            zodiac_deg_dms: moon_zodiac_dms,
            horiz_parallax_arcsec: dm.horiz_parallax_rad * RTS,
            semidiameter_arcsec: dm.semidiameter_arcsec,
            visual_magnitude: moon_mag,
            illuminated_fraction: dm.illuminated_fraction,
            solar_elongation_deg: dm.solar_elongation_deg,
            phase_note: dm.phase_description.clone(),
            topo_ra_hms: dm.topo_ra_hms.clone(),
            topo_dec_dms: dm.topo_dec_dms.clone(),
            topo_alt_deg: dm.topo_alt_deg,
            topo_az_deg: dm.topo_az_deg,
            refracted_alt_deg: dm.refracted_alt_deg,
            lha_hours: dm.lha_hours,
            rise_transit_set: dm.rise_transit_set.clone(),
        },
    );

    // Lunar Nodes & Lilith (`mean_elements` + osculating)
    let (mean_lilith, true_lilith, mean_node, true_node) = compute_lunar_nodes_and_lilith(jd_tdt);
    bodies.push(make_ecliptic_point_body(
        11,
        "Mean North Node (Rahu)",
        "gplan.c (mean_elements #17)",
        mean_node,
        jd_tdt,
        jd_ut,
        input.obs_lon_deg,
        input.obs_lat_deg,
        input.obs_height_m,
        input.temp_c,
        input.pressure_mb,
    ));
    bodies.push(make_ecliptic_point_body(
        12,
        "True North Node (Rahu)",
        "Osculating Lunar Node",
        true_node,
        jd_tdt,
        jd_ut,
        input.obs_lon_deg,
        input.obs_lat_deg,
        input.obs_height_m,
        input.temp_c,
        input.pressure_mb,
    ));
    bodies.push(make_ecliptic_point_body(
        13,
        "South Node (Ketu)",
        "Osculating Lunar South Node",
        modtp(true_node + std::f64::consts::PI),
        jd_tdt,
        jd_ut,
        input.obs_lon_deg,
        input.obs_lat_deg,
        input.obs_height_m,
        input.temp_c,
        input.pressure_mb,
    ));
    bodies.push(make_ecliptic_point_body(
        14,
        "Mean Lilith (Lunar Apogee)",
        "gplan.c (L - l + 180°)",
        modtp(mean_lilith + std::f64::consts::PI),
        jd_tdt,
        jd_ut,
        input.obs_lon_deg,
        input.obs_lat_deg,
        input.obs_height_m,
        input.temp_c,
        input.pressure_mb,
    ));
    bodies.push(make_ecliptic_point_body(
        15,
        "True Black Moon Lilith",
        "Osculating Lunar Apogee",
        modtp(true_lilith + std::f64::consts::PI),
        jd_tdt,
        jd_ut,
        input.obs_lon_deg,
        input.obs_lat_deg,
        input.obs_height_m,
        input.temp_c,
        input.pressure_mb,
    ));

    // `conjunct.c` summary (2 synodic cycles + 4 solar seasons)
    let conjunct_summary = compute_conjunct_report(jd_tdt, 2);

    // `moonrise.c` 10-day preview
    let moonrise_preview = compute_moonrise_table(
        jd_ut,
        10,
        input.obs_lon_deg,
        input.obs_lat_deg,
        input.obs_height_m,
        input.temp_c,
        input.pressure_mb,
    );

    EphemerisReport {
        engine: "aa56-rust (Moshier Astronomical Almanac v5.6 DE404 Native)".to_string(),
        jd_ut,
        jd_tdt,
        delta_t_sec,
        ut_date_str: jtocal_string(jd_ut),
        tdt_date_str: jtocal_string(jd_tdt),
        gst_hours: gst_rad * 12.0 / std::f64::consts::PI,
        gst_hms: format_hms_rad(gst_rad),
        lst_hours: lst_rad * 12.0 / std::f64::consts::PI,
        lst_hms: format_hms_rad(lst_rad),
        mean_obliquity_deg: eps * RTD,
        mean_obliquity_dms: format_dms(eps),
        nutation_lon_arcsec: nutl * RTS,
        nutation_obliq_arcsec: nuto * RTS,
        obs_lon_deg: input.obs_lon_deg,
        obs_lat_deg: input.obs_lat_deg,
        obs_height_m: input.obs_height_m,
        temp_c: input.temp_c,
        pressure_mb: input.pressure_mb,
        earth_helio_j2000,
        earth_helio_lon_deg: earth_polar_j2000[0] * RTD,
        earth_helio_lat_deg: earth_polar_j2000[1] * RTD,
        earth_helio_dist_au: earth_polar_j2000[2],
        bodies,
        domoon_details: dm,
        conjunct_summary,
        moonrise_preview,
    }
}
