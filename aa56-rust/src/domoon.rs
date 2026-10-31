//! Geocentric & Topocentric Moon Ephemeris (`domoon.c`).
//! Implements `moonll` and `domoon` with full DE404 lunar theory (`mlr404`, `mlat404`),
//! light-time aberration, horizontal parallax, semidiameter, phase angle, illuminated fraction,
//! quarter phase age, and topocentric reductions.

use serde::Serialize;

use crate::earth_orientation::{
    epsiln, geodetic_to_geocentric, lonlat, nutate, nutlo, precess, topocentric_altaz,
};
use crate::gplan::gmoon;
use crate::kepler::{kepler, De404Body, OrbitElement};
use crate::math_time::{
    deltap, format_dms, format_hms_rad, modtp, zatan2, AEARTH_M, AU_KM, DTR, PI, RTD, RTS, TPI,
};
use crate::transit::{compute_rise_transit_set, RiseSetTransit, HORIZON_MOON};

/// Earth equatorial radius in AU (`0.001 * aearth / au` in `domoon.c`).
pub const REARTH_AU: f64 = 0.001 * AEARTH_M / AU_KM;

#[derive(Debug, Clone, Serialize)]
pub struct DomoonResult {
    pub geo_rect_of_date: [f64; 3],
    pub geo_polar_of_date: [f64; 3],
    pub distance_earth_radii: f64,
    pub distance_au: f64,
    pub distance_km: f64,
    pub horiz_parallax_rad: f64,
    pub horiz_parallax_deg: f64,
    pub horiz_parallax_dms: String,
    pub semidiameter_rad: f64,
    pub semidiameter_arcsec: f64,
    pub semidiameter_dms: String,
    pub app_ra_rad: f64,
    pub app_dec_rad: f64,
    pub app_ra_hms: String,
    pub app_dec_dms: String,
    pub app_ecl_lon_deg: f64,
    pub app_ecl_lat_deg: f64,
    pub app_ecl_lon_dms: String,
    pub app_ecl_lat_dms: String,
    pub nutation_lon_arcsec: f64,
    pub nutation_obliq_arcsec: f64,
    pub solar_elongation_deg: f64,
    pub illuminated_fraction: f64,
    pub moon_sun_lon_diff_deg: f64,
    pub nearest_quarter_name: String,
    pub days_from_quarter: f64,
    pub phase_description: String,
    pub topo_ra_rad: f64,
    pub topo_dec_rad: f64,
    pub topo_ra_hms: String,
    pub topo_dec_dms: String,
    pub topo_alt_deg: f64,
    pub topo_az_deg: f64,
    pub refracted_alt_deg: f64,
    pub lha_hours: f64,
    pub rise_transit_set: RiseSetTransit,
}

/// Calculate apparent latitude, longitude, and horizontal parallax of the Moon at Julian date `jed` (TDT)
/// (`moonll` in `domoon.c`).
/// Returns `(rect_apparent_of_date, pol_apparent_of_date, angles_ep, angles_pq, angles_ dlon_rad)`.
pub fn moonll(jed: f64, rearth_j2000: &[f64; 3]) -> ([f64; 3], [f64; 3], f64, f64, f64) {
    let eps = epsiln(jed);
    let coseps = eps.cos();
    let sineps = eps.sin();

    // Compute geocentric Moon coordinates (equatorial rect of date, ecliptic polar of date)
    let (_rect0, mut pol) = gmoon(jed);

    // Light-time correction to longitude (`0.0118 * DTR * Rearth / pol[2]` in `domoon.c`)
    pol[0] -= 0.0118 * DTR * REARTH_AU / pol[2];

    // Convert to equatorial rectangular of date
    let cosB = pol[1].cos();
    let sinB = pol[1].sin();
    let cosL = pol[0].cos();
    let sinL = pol[0].sin();
    let x_ec = cosB * cosL * pol[2];
    let y_ec = cosB * sinL * pol[2];
    let z_ec = sinB * pol[2];

    let mut rect = [
        x_ec,
        coseps * y_ec - sineps * z_ec,
        sineps * y_ec + coseps * z_ec,
    ];

    // Rotate to J2000 (`precess(rect, J, 1)`)
    precess(&mut rect, jed, 1);

    // Find Euclidean vectors and angles between earth, object, and the sun (`angles` in `domoon.c`)
    let mut pp = [0.0f64; 3];
    let mut qq = [0.0f64; 3];
    for i in 0..3 {
        pp[i] = rect[i];
        qq[i] = rearth_j2000[i] + pp[i];
    }
    let (ep, pq) = compute_angles(&pp, &qq, rearth_j2000);

    // Make rect a unit vector
    let mut r_norm = 0.0;
    for i in 0..3 {
        r_norm += rect[i] * rect[i];
    }
    r_norm = r_norm.sqrt();
    for i in 0..3 {
        rect[i] /= r_norm;
    }

    // Correct position for light deflection? (`domoon.c` leaves `relativity` commented out)
    // Precess to equinox of date
    precess(&mut rect, jed, -1);

    // Geocentric ecliptic longitude and latitude of date (`lonlat` in `domoon.c`)
    let r_dist = pol[2];
    let (lon_rad, lat_rad, _) = lonlat(&rect, jed);
    pol[0] = lon_rad;
    pol[1] = lat_rad;
    pol[2] = r_dist;

    // Nutation (`nutate(J, rect)`)
    let (_dpsi, _deps) = nutate(jed, &mut rect);

    // Sun ecliptic longitude of date for phase quarter calculation
    let mut sun_neg = [-rearth_j2000[0], -rearth_j2000[1], -rearth_j2000[2]];
    precess(&mut sun_neg, jed, -1);
    let (sun_lon, _, _) = lonlat(&sun_neg, jed);
    let dlon = deltap(pol[0], sun_lon);

    (rect, pol, ep, pq, dlon)
}

/// `angles(p, q, e)` from `angles.c` returning `(ep, pq)` cosine angles.
fn compute_angles(p: &[f64; 3], q: &[f64; 3], e: &[f64; 3]) -> (f64, f64) {
    let mut a = 0.0;
    let mut b = 0.0;
    let mut c = 0.0;
    let mut ep = 0.0;
    let mut pq = 0.0;
    for i in 0..3 {
        a += e[i] * e[i];
        b += p[i] * p[i];
        c += q[i] * q[i];
        ep += e[i] * p[i];
        pq += p[i] * q[i];
    }
    let eo = a.sqrt();
    let po = b.sqrt();
    let qo = c.sqrt();
    let ep_cos = if eo * po > 0.0 { ep / (eo * po) } else { 0.0 };
    let pq_cos = if po * qo > 0.0 { pq / (po * qo) } else { 0.0 };
    (ep_cos, pq_cos)
}

/// Complete `domoon` calculation (`domoon.c`).
pub fn domoon(
    jd_tdt: f64,
    jd_ut: f64,
    obs_lon_deg: f64,
    obs_lat_deg: f64,
    obs_height_m: f64,
    temp_c: f64,
    pressure_mb: f64,
) -> DomoonResult {
    // 1. Ensure Earth heliocentric J2000 position is computed at TDT
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
    let (rearth_j2000, _) = kepler(jd_tdt, &mut earth_el);

    // 2. Run `moonll(TDT, moonpp, moonpol)`
    let (moonpp, moonpol, ep, pq, dlon) = moonll(jd_tdt, &rearth_j2000);
    let (nutl, nuto) = nutlo(jd_tdt);

    // Apparent RA & Dec (`zatan2(moonpp[0], moonpp[1])`, `asin(moonpp[2])`)
    let ra = zatan2(moonpp[0], moonpp[1]);
    let dec = moonpp[2].clamp(-1.0, 1.0).asin();

    // Horizontal parallax (`asin(Rearth / moonpol[2])`)
    let hp_sin = (REARTH_AU / moonpol[2]).clamp(-1.0, 1.0);
    let obj_hp = hp_sin.asin();

    // Semidiameter (`0.272453 * x + 0.0799 / RTS` in `domoon.c`)
    let obj_sd = 0.272453 * hp_sin + 0.0799 / RTS;

    let dist_earth_radii = moonpol[2] / REARTH_AU;
    let dist_au = moonpol[2];
    let dist_km = dist_au * AU_KM;

    // Elongation from Sun (`RTD * acos(-ep)`) and Illuminated fraction (`0.5 * (1.0 + pq)`)
    let solar_elongation_deg = RTD * (-ep).clamp(-1.0, 1.0).acos();
    let illuminated_fraction = 0.5 * (1.0 + pq);

    // Quarter phase & age (`domoon.c` lines 163-172)
    let x_phase = dlon / (0.5 * PI);
    let i_quarter = if x_phase >= 0.0 {
        (x_phase + 0.5) as i32
    } else {
        (x_phase - 0.5) as i32
    };
    let days_from_quarter = (dlon - (i_quarter as f64) * 0.5 * PI) / (DTR * 12.19075);
    let q_idx = ((i_quarter + 2) & 3) as usize;
    let quarters = ["Full Moon", "Third Quarter", "New Moon", "First Quarter"];
    let nearest_quarter_name = quarters[q_idx].to_string();
    let phase_description = format!(
        "Moon phase {:.1} days {} {}",
        days_from_quarter.abs(),
        if days_from_quarter < 0.0 {
            "before"
        } else {
            "past"
        },
        nearest_quarter_name
    );

    // Topocentric reduction (`altaz.c`)
    let geo_obs = geodetic_to_geocentric(obs_lat_deg, obs_height_m);
    let topo = topocentric_altaz(
        ra,
        dec,
        jd_ut,
        jd_tdt,
        obs_lon_deg,
        obs_lat_deg,
        &geo_obs,
        obj_hp,
        temp_c,
        pressure_mb,
    );

    // Iterative Rise, Transit, Set (`iter_trnsit` in `domoon.c`)
    let dt_days = jd_tdt - jd_ut;
    let rts = compute_rise_transit_set(
        jd_ut,
        obs_lon_deg,
        obs_lat_deg,
        &geo_obs,
        HORIZON_MOON,
        temp_c,
        pressure_mb,
        |t_ut| {
            let t_tdt = t_ut + dt_days;
            let (re, _) = kepler(t_tdt, &mut earth_el);
            let (mpp, mpol, _, _, _) = moonll(t_tdt, &re);
            let r_ra = zatan2(mpp[0], mpp[1]);
            let r_dec = mpp[2].clamp(-1.0, 1.0).asin();
            let r_hp = (REARTH_AU / mpol[2]).clamp(-1.0, 1.0).asin();
            let r_sd = 0.272453 * (REARTH_AU / mpol[2]) + 0.0799 / RTS;
            (r_ra, r_dec, r_hp, r_sd)
        },
    );

    DomoonResult {
        geo_rect_of_date: moonpp,
        geo_polar_of_date: moonpol,
        distance_earth_radii: dist_earth_radii,
        distance_au: dist_au,
        distance_km: dist_km,
        horiz_parallax_rad: obj_hp,
        horiz_parallax_deg: obj_hp * RTD,
        horiz_parallax_dms: format_dms(obj_hp),
        semidiameter_rad: obj_sd,
        semidiameter_arcsec: obj_sd * RTS,
        semidiameter_dms: format_dms(obj_sd),
        app_ra_rad: ra,
        app_dec_rad: dec,
        app_ra_hms: format_hms_rad(ra),
        app_dec_dms: format_dms(dec),
        app_ecl_lon_deg: modtp(moonpol[0]) * RTD,
        app_ecl_lat_deg: moonpol[1] * RTD,
        app_ecl_lon_dms: format_dms(modtp(moonpol[0])),
        app_ecl_lat_dms: format_dms(moonpol[1]),
        nutation_lon_arcsec: nutl * RTS,
        nutation_obliq_arcsec: nuto * RTS,
        solar_elongation_deg,
        illuminated_fraction,
        moon_sun_lon_diff_deg: dlon * RTD,
        nearest_quarter_name,
        days_from_quarter,
        phase_description,
        topo_ra_rad: topo.topo_ra_rad,
        topo_dec_rad: topo.topo_dec_rad,
        topo_ra_hms: format_hms_rad(topo.topo_ra_rad),
        topo_dec_dms: format_dms(topo.topo_dec_rad),
        topo_alt_deg: topo.true_alt_deg,
        topo_az_deg: topo.true_az_deg,
        refracted_alt_deg: topo.refracted_alt_deg,
        lha_hours: topo.lha_hours,
        rise_transit_set: rts,
    }
}

