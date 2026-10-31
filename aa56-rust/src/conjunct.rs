//! Search for exterior or interior conjunctions of Moon and Sun, Lunar Quarters, and Equinoxes/Solstices (`conjunct.c`).
//! Implements `apparent`, `func`, and `search` from Stephen L. Moshier's `conjunct.c`.

use serde::Serialize;

use crate::domoon::moonll;
use crate::earth_orientation::{epsiln, lonlat, nutate, precess};
use crate::kepler::{kepler, De404Body, OrbitElement};
use crate::math_time::{deltap, deltat, jtocal_string, zatan2, J2000, PI, RTD, STR, TPI};
use crate::reductions::annuab;

#[derive(Debug, Clone, Serialize)]
pub struct ConjunctionEvent {
    pub event_type: String,
    pub event_category: String,
    pub target_angle_deg: f64,
    pub jd_tdt: f64,
    pub jd_ut: f64,
    pub date_tdt_str: String,
    pub date_ut_str: String,
    pub sun_app_lon_deg: f64,
    pub moon_app_lon_deg: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct ConjunctReport {
    pub center_jd: f64,
    pub lunar_phases: Vec<ConjunctionEvent>,
    pub solar_seasons: Vec<ConjunctionEvent>,
}

/// Compute apparent ecliptic longitudes of the Sun and Moon at Julian date `jd_tdt` (`apparent` in `conjunct.c`).
/// Returns `(sun_app_lon_rad, moon_app_lon_rad)` in `[0, 2*PI)`.
pub fn apparent_sun_moon_lon(jd_tdt: f64) -> (f64, f64) {
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

    // Present position of Earth
    let (rearth, eapolar) = kepler(jd_tdt, &mut earth_el);

    // Apparent position of the Moon (`moonll`)
    let (_moonpp, moonpol, _, _, _) = moonll(jd_tdt, &rearth);
    let m = moonpol[0];

    // Apparent position of the Sun: include light time (`t - 0.0057755183 * eapolar[2]`)
    let r_earth = eapolar[2];
    let (ecr_past, _) = kepler(jd_tdt - 0.0057755183 * r_earth, &mut earth_el);
    let mut s_vec = [-ecr_past[0], -ecr_past[1], -ecr_past[2]];

    // Annual aberration (`annuab`)
    annuab(&mut s_vec, jd_tdt);
    // Precess from J2000 to date (`precess(s, t, -1)`)
    precess(&mut s_vec, jd_tdt, -1);
    // Ecliptic longitude of date (`lonlat(s, t, polar, 0)`)
    let (s_lon, _, _) = lonlat(&s_vec, jd_tdt);

    (s_lon, m)
}

/// Search for date when `func(t) == y` using `search(t, y, delta)` from `conjunct.c`.
pub fn search_conjunct<F>(t_guess: f64, y_target: f64, delta: f64, func: F) -> f64
where
    F: Fn(f64) -> f64,
{
    let mut tl = t_guess - 0.5 * delta;
    let mut tm = t_guess;
    let mut tr = t_guess + 0.5 * delta;

    let mut yl = deltap(func(tl), y_target);
    let mut ym = deltap(func(tm), y_target);
    let mut yr = deltap(func(tr), y_target);

    for _ in 0..16 {
        let a = 0.5 * (yl + yr) - ym;
        let b = 0.5 * (yr - yl);
        let c = ym;
        let d = b * b - 4.0 * a * c;

        let e = if d < 0.0 || a.abs() < 1.0e-14 {
            if b.abs() > 1.0e-14 {
                -c / b
            } else {
                0.0
            }
        } else if b < 0.0 {
            (-b - d.sqrt()) / (2.0 * a)
        } else {
            (-b + d.sqrt()) / (2.0 * a)
        };

        let tn = tm + 0.5 * e * (tr - tl);
        let yn = deltap(func(tn), y_target);
        if yn.abs() < STR * 0.01 || (tr - tl).abs() < 1.0e-8 {
            return tn;
        }

        if e < 0.0 {
            tr = tm;
            yr = ym;
        } else {
            tl = tm;
            yl = ym;
        }
        tm = tn;
        ym = yn;
    }
    tm
}

/// Generate a complete report of Lunar Phases (New, First Quarter, Full, Third Quarter)
/// and Solar Equinoxes/Solstices around `jd_start` (`conjunct.c`).
pub fn compute_conjunct_report(jd_start: f64, num_synodic_cycles: usize) -> ConjunctReport {
    let mut lunar_phases = Vec::new();
    let synodic_month = 29.530588853;

    // Find current elongation `s - m` at `jd_start`
    let (s0, m0) = apparent_sun_moon_lon(jd_start);
    let diff0 = deltap(m0, s0); // `m - s` in `[-PI, PI]`
    let mut approx_new_moon = jd_start - (diff0 / TPI) * synodic_month;

    // Step back 1 cycle if ahead of `jd_start` so we cover starting from just before `jd_start`
    if approx_new_moon > jd_start {
        approx_new_moon -= synodic_month;
    }

    // In `conjunct.c`, `func(t)` returns `s - m`:
    // New Moon: `s - m == 0.0`
    // First Quarter: `s - m == -0.5 * PI` (`m - s == +90°`)
    // Full Moon: `s - m == PI` (`m - s == 180°`)
    // Third Quarter: `s - m == 0.5 * PI` (`m - s == 270°`)
    let phase_defs = [
        ("New Moon (Lua Nova)", 0.0, 0.0, 0.0),
        ("First Quarter (Quarto Crescente)", -0.5 * PI, 90.0, 0.25),
        ("Full Moon (Lua Cheia)", PI, 180.0, 0.50),
        ("Third Quarter (Quarto Minguante)", 0.5 * PI, 270.0, 0.75),
    ];

    for cycle in 0..num_synodic_cycles.max(1) {
        let base_t = approx_new_moon + (cycle as f64) * synodic_month;
        for &(name, target_s_minus_m, deg_label, offset_frac) in &phase_defs {
            let guess = base_t + offset_frac * synodic_month;
            let t_tdt = search_conjunct(guess, target_s_minus_m, 2.0, |t| {
                let (s, m) = apparent_sun_moon_lon(t);
                s - m
            });
            let (s_rad, m_rad) = apparent_sun_moon_lon(t_tdt);
            let y_approx = 2000.0 + (t_tdt - J2000) / 365.2425;
            let dt_sec = deltat(y_approx);
            let t_ut = t_tdt - dt_sec / 86400.0;

            lunar_phases.push(ConjunctionEvent {
                event_type: name.to_string(),
                event_category: "Lunar Phase".to_string(),
                target_angle_deg: deg_label,
                jd_tdt: t_tdt,
                jd_ut: t_ut,
                date_tdt_str: jtocal_string(t_tdt),
                date_ut_str: jtocal_string(t_ut),
                sun_app_lon_deg: s_rad * RTD,
                moon_app_lon_deg: m_rad * RTD,
            });
        }
    }

    // Also compute the 4 Solar Seasons (Vernal Equinox 0°, Summer Solstice 90°, Autumn Equinox 180°, Winter Solstice 270°)
    // around the year of `jd_start` (`conjunct.c` #if EQUINOX).
    let mut solar_seasons = Vec::new();
    let season_defs = [
        ("March Equinox (Equinócio λ☉ = 0°)", 0.0, 0.0),
        ("June Solstice (Solstício λ☉ = 90°)", 0.5 * PI, 90.0),
        ("September Equinox (Equinócio λ☉ = 180°)", PI, 180.0),
        ("December Solstice (Solstício λ☉ = 270°)", 1.5 * PI, 270.0),
    ];

    // Estimate March equinox near the year of `jd_start`
    let s0_pos = if s0 < 0.0 { s0 + TPI } else { s0 };
    let approx_vernal = jd_start - (s0_pos / TPI) * 365.2422;

    for &(name, target_s, deg_label) in &season_defs {
        let guess = approx_vernal + (target_s / TPI) * 365.2422;
        let t_tdt = search_conjunct(guess, target_s, 3.0, |t| {
            let (s, _m) = apparent_sun_moon_lon(t);
            s
        });
        let (s_rad, m_rad) = apparent_sun_moon_lon(t_tdt);
        let y_approx = 2000.0 + (t_tdt - J2000) / 365.2425;
        let dt_sec = deltat(y_approx);
        let t_ut = t_tdt - dt_sec / 86400.0;

        solar_seasons.push(ConjunctionEvent {
            event_type: name.to_string(),
            event_category: "Solar Season".to_string(),
            target_angle_deg: deg_label,
            jd_tdt: t_tdt,
            jd_ut: t_ut,
            date_tdt_str: jtocal_string(t_tdt),
            date_ut_str: jtocal_string(t_ut),
            sun_app_lon_deg: s_rad * RTD,
            moon_app_lon_deg: m_rad * RTD,
        });
    }

    ConjunctReport {
        center_jd: jd_start,
        lunar_phases,
        solar_seasons,
    }
}
