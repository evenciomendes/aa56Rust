//! Star, Messier, and Orbit Catalogs Parser & Reductions (`rstar.c`, `fk4fk5.c`, `kfiles.c`).
//! Supports `star.cat`, `star1950.cat`, `messier.cat`, and `orbit.cat` from `aa-56`.

use serde::Serialize;

use crate::constel::{expand_star_name, whatconstel};
use crate::earth_orientation::{
    geodetic_to_geocentric, lonlat, nutate, nutlo, precess, topocentric_altaz,
};
use crate::kepler::{kepler, De404Body, OrbitElement};
use crate::math_time::{
    caltoj, format_dms, format_hms_rad, modtp, zatan2, B1950, CLIGHT_KMS, DTR, J2000, RTD, RTS,
    STR,
};
use crate::reductions::{annuab, format_zodiac, reduce_planet, relativity, BodyEphemeris};
use crate::transit::{compute_rise_transit_set, RiseSetTransit, HORIZON_PLANET_STAR};

#[derive(Debug, Clone)]
pub struct StarCatalogRecord {
    pub raw_name: String,
    pub expanded_name: String,
    pub epoch: f64,
    pub ra_rad: f64,
    pub dec_rad: f64,
    pub mu_ra_sec_century: f64,
    pub mu_dec_arcsec_century: f64,
    pub v_rad_kms: f64,
    pub parallax_arcsec: f64,
    pub visual_mag: f64,
    pub description: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct ReducedStarEntry {
    pub index: usize,
    pub raw_name: String,
    pub expanded_name: String,
    pub catalog_epoch: f64,
    pub constellation: String,
    pub visual_mag: f64,
    pub parallax_arcsec: f64,
    pub radial_vel_kms: f64,
    pub proper_motion_ra_s_cy: f64,
    pub proper_motion_dec_as_cy: f64,
    pub j2000_ra_hms: String,
    pub j2000_dec_dms: String,
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
    pub topo_ra_hms: String,
    pub topo_dec_dms: String,
    pub topo_alt_deg: f64,
    pub topo_az_deg: f64,
    pub refracted_alt_deg: f64,
    pub lha_hours: f64,
    pub rise_transit_set: RiseSetTransit,
    pub description: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct ReducedOrbitEntry {
    pub index: usize,
    pub epoch: f64,
    pub inclination_deg: f64,
    pub asc_node_deg: f64,
    pub arg_perihelion_deg: f64,
    pub semi_major_axis_au: f64,
    pub daily_motion_deg: f64,
    pub eccentricity: f64,
    pub mean_anomaly_deg: f64,
    pub equinox: f64,
    pub abs_mag: f64,
    pub ephemeris: BodyEphemeris,
}

/// FK4 (B1950) to FK5 (J2000) transformation (`fk4fk5.c` by Yallop et al. 1989).
pub fn fk4fk5(star: &mut StarCatalogRecord) {
    let m_mat: [[f64; 6]; 6] = [
        [
            0.9999256782,
            -0.0111820611,
            -0.0048579477,
            0.00000242395018,
            -0.00000002710663,
            -0.00000001177656,
        ],
        [
            0.0111820610,
            0.9999374784,
            -0.0000271765,
            0.00000002710663,
            0.00000242397878,
            -0.00000000006587,
        ],
        [
            0.0048579479,
            -0.0000271474,
            0.9999881997,
            0.00000001177656,
            -0.00000000006582,
            0.00000242410173,
        ],
        [
            -0.000551,
            -0.238565,
            0.435739,
            0.99994704,
            -0.01118251,
            -0.00485767,
        ],
        [
            0.238514,
            -0.002667,
            -0.008541,
            0.01118251,
            0.99995883,
            -0.00002718,
        ],
        [
            -0.435623,
            0.012254,
            0.002117,
            0.00485767,
            -0.00002714,
            1.00000956,
        ],
    ];

    let a_vec: [f64; 6] = [
        -1.62557e-6,
        -0.31919e-6,
        -0.13843e-6,
        1.245e-3,
        -1.580e-3,
        -0.659e-3,
    ];

    let a = star.ra_rad;
    let d = star.dec_rad;
    let sina = a.sin();
    let cosa = a.cos();
    let sind = d.sin();
    let cosd = d.cos();

    let r0 = [cosd * cosa, cosd * sina, sind];

    // Proper motion in FK4 is in seconds of time/century (RA) and arcsec/century (Dec)
    let mura = 15.0 * star.mu_ra_sec_century;
    let mudec = star.mu_dec_arcsec_century;
    let px = star.parallax_arcsec;
    let v = 21.095 * star.v_rad_kms * px;

    let r1 = [
        -mura * cosd * sina - mudec * sind * cosa + v * r0[0],
        mura * cosd * cosa - mudec * sind * sina + v * r0[1],
        mudec * cosd + v * r0[2],
    ];

    let mut pv0 = [0.0f64; 6];
    for i in 0..3 {
        pv0[i] = r0[i] - a_vec[i];
        pv0[i + 3] = r1[i] - a_vec[i + 3];
    }

    let mut pv1 = [0.0f64; 6];
    for i in 0..6 {
        let mut s = 0.0;
        for j in 0..6 {
            s += m_mat[i][j] * pv0[j];
        }
        pv1[i] = s;
    }

    let x = pv1[0];
    let y = pv1[1];
    let z = pv1[2];
    let r = (x * x + y * y + z * z).sqrt();
    let xd = pv1[3] / r;
    let yd = pv1[4] / r;
    let zd = pv1[5] / r;
    let xp = x / r;
    let yp = y / r;
    let zp = z / r;

    let ra_new = zatan2(xp, yp);
    let dec_new = zp.clamp(-1.0, 1.0).asin();

    let rxy2 = xp * xp + yp * yp;
    let rd = xp * xd + yp * yd + zp * zd;

    let mura_new = if rxy2 > 1.0e-14 {
        (xp * yd - yp * xd) / rxy2
    } else {
        0.0
    };
    let mudec_new = if rxy2 > 1.0e-14 {
        (zd - zp * rd) / rxy2.sqrt()
    } else {
        0.0
    };

    if px > 0.0 {
        star.v_rad_kms = rd / (21.095 * px);
        star.parallax_arcsec = px / r;
    }
    star.ra_rad = ra_new;
    star.dec_rad = dec_new;
    star.mu_ra_sec_century = mura_new / 15.0;
    star.mu_dec_arcsec_century = mudec_new;
    star.epoch = J2000;
}

/// Parse a star or Messier catalog file (`kfiles.c` -> `getstar`).
pub fn parse_star_catalog(content: &str) -> Vec<StarCatalogRecord> {
    let mut records = Vec::new();
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with(';') {
            continue;
        }
        let tokens: Vec<&str> = trimmed.split_whitespace().collect();
        if tokens.len() < 12 {
            continue;
        }

        // Format 1 (`star.cat`, `star1950.cat`):
        // epoch rh rm rs dd dm ds mura mudec v px mag name [desc...]
        // Format 2 (`messier.cat`):
        // rh rm rs dd dm ds mura mudec v px mag epoch name [desc...]
        let first_val = tokens[0].parse::<f64>().unwrap_or(0.0);
        let (
            epoch_raw,
            rh,
            rm,
            rs,
            dd_str,
            dm,
            ds,
            mura,
            mudec,
            v,
            px,
            mag,
            name_tok,
            desc_str,
        ) = if first_val > 1000.0 {
            let epoch = first_val;
            let rh = tokens[1].parse::<f64>().unwrap_or(0.0);
            let rm = tokens[2].parse::<f64>().unwrap_or(0.0);
            let rs = tokens[3].parse::<f64>().unwrap_or(0.0);
            let dd_str = tokens[4];
            let dm = tokens[5].parse::<f64>().unwrap_or(0.0);
            let ds = tokens[6].parse::<f64>().unwrap_or(0.0);
            let mura = tokens[7].parse::<f64>().unwrap_or(0.0);
            let mudec = tokens[8].parse::<f64>().unwrap_or(0.0);
            let v = tokens[9].parse::<f64>().unwrap_or(0.0);
            let px = tokens[10].parse::<f64>().unwrap_or(0.0);
            let mag = tokens[11].parse::<f64>().unwrap_or(0.0);
            let name_tok = if tokens.len() > 12 { tokens[12] } else { "Object" };
            let desc = if tokens.len() > 13 {
                tokens[13..].join(" ")
            } else {
                String::new()
            };
            (epoch, rh, rm, rs, dd_str, dm, ds, mura, mudec, v, px, mag, name_tok, desc)
        } else {
            let rh = first_val;
            let rm = tokens[1].parse::<f64>().unwrap_or(0.0);
            let rs = tokens[2].parse::<f64>().unwrap_or(0.0);
            let dd_str = tokens[3];
            let dm = tokens[4].parse::<f64>().unwrap_or(0.0);
            let ds = tokens[5].parse::<f64>().unwrap_or(0.0);
            let mura = tokens[6].parse::<f64>().unwrap_or(0.0);
            let mudec = tokens[7].parse::<f64>().unwrap_or(0.0);
            let v = tokens[8].parse::<f64>().unwrap_or(0.0);
            let px = tokens[9].parse::<f64>().unwrap_or(0.0);
            let mag = tokens[10].parse::<f64>().unwrap_or(0.0);
            let epoch = tokens[11].parse::<f64>().unwrap_or(2000.0);
            let name_tok = if tokens.len() > 12 { tokens[12] } else { "Messier" };
            let desc = if tokens.len() > 13 {
                tokens[13..].join(" ")
            } else {
                String::new()
            };
            (epoch, rh, rm, rs, dd_str, dm, ds, mura, mudec, v, px, mag, name_tok, desc)
        };

        let epoch_jd = if (epoch_raw - 2000.0).abs() < 1e-3 {
            J2000
        } else if (epoch_raw - 1950.0).abs() < 1e-3 {
            B1950
        } else if epoch_raw < 10000.0 {
            caltoj(epoch_raw as i64, 1, 1.5)
        } else {
            epoch_raw
        };

        let ra_rad = 15.0 * DTR * (rh + rm / 60.0 + rs / 3600.0);
        let dec_sign = if dd_str.starts_with('-') { -1.0 } else { 1.0 };
        let dd_abs = dd_str.parse::<f64>().unwrap_or(0.0).abs();
        let dec_rad = dec_sign * DTR * (dd_abs + dm / 60.0 + ds / 3600.0);

        let px_arcsec = if px >= 1.0 { 0.001 * px } else { px };
        let expanded_name = expand_star_name(name_tok);

        records.push(StarCatalogRecord {
            raw_name: name_tok.to_string(),
            expanded_name,
            epoch: epoch_jd,
            ra_rad,
            dec_rad,
            mu_ra_sec_century: mura,
            mu_dec_arcsec_century: mudec,
            v_rad_kms: v,
            parallax_arcsec: px_arcsec,
            visual_mag: mag,
            description: desc_str,
        });
    }
    records
}

/// Reduce a fixed star or Messier object (`rstar.c` -> `doel`).
pub fn reduce_star(
    index: usize,
    rec: &StarCatalogRecord,
    earth_j2000: &[f64; 3],
    jd_tdt: f64,
    jd_ut: f64,
    obs_lon_deg: f64,
    obs_lat_deg: f64,
    obs_height_m: f64,
    temp_c: f64,
    pressure_mb: f64,
) -> ReducedStarEntry {
    let mut star = rec.clone();
    if (star.epoch - B1950).abs() < 1.0 {
        fk4fk5(&mut star);
    }

    let cosdec = star.dec_rad.cos();
    let sindec = star.dec_rad.sin();
    let cosra = star.ra_rad.cos();
    let sinra = star.ra_rad.sin();

    let mut q = [cosdec * cosra, cosdec * sinra, sindec];

    // Space motion (proper motion + radial velocity + parallax in `rstar.c`)
    let epoch_diff_days = jd_tdt - star.epoch;
    let mura_rad_day = (15.0 * star.mu_ra_sec_century * STR) / 36525.0;
    let mudec_rad_day = (star.mu_dec_arcsec_century * STR) / 36525.0;
    let vpi_day = 21.094952663 * star.v_rad_kms * star.parallax_arcsec / (36525.0 * RTS);

    let m = [
        (-mura_rad_day * cosdec * sinra - mudec_rad_day * sindec * cosra + vpi_day * q[0])
            * epoch_diff_days,
        (mura_rad_day * cosdec * cosra - mudec_rad_day * sindec * sinra + vpi_day * q[1])
            * epoch_diff_days,
        (mudec_rad_day * cosdec + vpi_day * q[2]) * epoch_diff_days,
    ];
    for i in 0..3 {
        q[i] += m[i];
    }

    // Precess to J2000 if catalog epoch was neither B1950 nor J2000
    if (star.epoch - J2000).abs() > 1.0 {
        precess(&mut q, star.epoch, 1);
    }

    // Annual parallax (`p[i] = q[i] - px * e[i]`)
    let px_rad = star.parallax_arcsec * STR;
    let mut p = [0.0f64; 3];
    let mut r = 0.0;
    for i in 0..3 {
        p[i] = q[i] - px_rad * earth_j2000[i];
        r += p[i] * p[i];
    }
    r = r.sqrt();
    for i in 0..3 {
        p[i] /= r;
    }

    let j2000_ra = zatan2(p[0], p[1]);
    let j2000_dec = p[2].clamp(-1.0, 1.0).asin();
    let constel = whatconstel(&p, J2000).to_string();

    // Relativistic light deflection (`relativity(p, q, e)` with `q` at 10000 AU)
    let q_far = [q[0] * 10000.0, q[1] * 10000.0, q[2] * 10000.0];
    relativity(&mut p, &q_far, earth_j2000);

    // Annual aberration (`annuab(p)`)
    annuab(&mut p, jd_tdt);

    // Precess from J2000 to date (`precess(p, TDT, -1)`)
    precess(&mut p, jd_tdt, -1);

    let (ecl_lon_mean, ecl_lat_mean, _) = lonlat(&p, jd_tdt);
    let (nutl, _) = nutlo(jd_tdt);
    let app_ecl_lon = modtp(ecl_lon_mean + nutl);
    let app_ecl_lat = ecl_lat_mean;

    // Nutation (`nutate(TDT, p)`)
    nutate(jd_tdt, &mut p);

    let app_ra = zatan2(p[0], p[1]);
    let app_dec = p[2].clamp(-1.0, 1.0).asin();
    let (zodiac_sign, zodiac_deg_dms) = format_zodiac(app_ecl_lon * RTD);

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

    ReducedStarEntry {
        index,
        raw_name: rec.raw_name.clone(),
        expanded_name: rec.expanded_name.clone(),
        catalog_epoch: if (rec.epoch - B1950).abs() < 1.0 {
            1950.0
        } else {
            2000.0
        },
        constellation: constel,
        visual_mag: rec.visual_mag,
        parallax_arcsec: rec.parallax_arcsec,
        radial_vel_kms: rec.v_rad_kms,
        proper_motion_ra_s_cy: rec.mu_ra_sec_century,
        proper_motion_dec_as_cy: rec.mu_dec_arcsec_century,
        j2000_ra_hms: format_hms_rad(j2000_ra),
        j2000_dec_dms: format_dms(j2000_dec),
        apparent_ra_rad: app_ra,
        apparent_dec_rad: app_dec,
        apparent_ra_hms: format_hms_rad(app_ra),
        apparent_dec_dms: format_dms(app_dec),
        apparent_ecl_lon_deg: app_ecl_lon * RTD,
        apparent_ecl_lat_deg: app_ecl_lat * RTD,
        apparent_ecl_lon_dms: format_dms(app_ecl_lon),
        apparent_ecl_lat_dms: format_dms(app_ecl_lat),
        zodiac_sign,
        zodiac_deg_dms,
        topo_ra_hms: format_hms_rad(topo.topo_ra_rad),
        topo_dec_dms: format_dms(topo.topo_dec_rad),
        topo_alt_deg: topo.true_alt_deg,
        topo_az_deg: topo.true_az_deg,
        refracted_alt_deg: topo.refracted_alt_deg,
        lha_hours: topo.lha_hours,
        rise_transit_set: rts,
        description: rec.description.clone(),
    }
}

/// Parse `orbit.cat` (`kfiles.c` -> `getorbit`).
pub fn parse_orbit_catalog(content: &str) -> Vec<OrbitElement> {
    let mut orbits = Vec::new();
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with(';') {
            continue;
        }
        let tokens: Vec<&str> = trimmed.split_whitespace().collect();
        if tokens.len() < 10 {
            continue;
        }
        let epoch_raw = tokens[0].parse::<f64>().unwrap_or(2451545.0);
        let i = tokens[1].parse::<f64>().unwrap_or(0.0);
        let cap_w = tokens[2].parse::<f64>().unwrap_or(0.0);
        let low_w = tokens[3].parse::<f64>().unwrap_or(0.0);
        let a = tokens[4].parse::<f64>().unwrap_or(1.0);
        let dm = tokens[5].parse::<f64>().unwrap_or(0.0);
        let ecc = tokens[6].parse::<f64>().unwrap_or(0.0);
        let m_anom = tokens[7].parse::<f64>().unwrap_or(0.0);
        let eq_raw = tokens[8].parse::<f64>().unwrap_or(2000.0);

        let (mag, sdiam, name) = if tokens.len() >= 11 {
            if let Ok(m_val) = tokens[9].parse::<f64>() {
                if let Ok(sd_val) = tokens[10].parse::<f64>() {
                    let nm = if tokens.len() > 11 {
                        tokens[11..].join(" ")
                    } else {
                        tokens[10].to_string()
                    };
                    (m_val, sd_val, nm)
                } else {
                    (m_val, 0.0, tokens[10..].join(" "))
                }
            } else {
                (8.0, 0.0, tokens[9..].join(" "))
            }
        } else {
            (8.0, 0.0, tokens[9].to_string())
        };

        let equinox = if (eq_raw - 2000.0).abs() < 1e-3 {
            J2000
        } else if (eq_raw - 1950.0).abs() < 1e-3 {
            B1950
        } else if eq_raw < 10000.0 {
            caltoj(eq_raw as i64, 1, 1.5)
        } else {
            eq_raw
        };

        let epoch = if epoch_raw < 10000.0 {
            caltoj(epoch_raw as i64, 1, 1.5)
        } else {
            epoch_raw
        };

        orbits.push(OrbitElement {
            obname: name,
            epoch,
            i,
            W: cap_w,
            w: low_w,
            a,
            dm,
            ecc,
            M: m_anom,
            equinox,
            mag,
            sdiam,
            L: 0.0,
            r: 0.0,
            plat: 0.0,
            de404: De404Body::None,
        });
    }
    orbits
}

/// Reduce all objects in `orbit.cat` at `jd_tdt`.
pub fn reduce_orbit_catalog(
    content: &str,
    earth_j2000: &[f64; 3],
    jd_tdt: f64,
    jd_ut: f64,
    obs_lon_deg: f64,
    obs_lat_deg: f64,
    obs_height_m: f64,
    temp_c: f64,
    pressure_mb: f64,
) -> Vec<ReducedOrbitEntry> {
    let mut elements = parse_orbit_catalog(content);
    let mut out = Vec::with_capacity(elements.len());
    for (idx, elem) in elements.iter_mut().enumerate() {
        let eph = reduce_planet(
            (idx + 100) as i32,
            "Asteroid / Comet (orbit.cat)",
            elem,
            earth_j2000,
            jd_tdt,
            jd_ut,
            obs_lon_deg,
            obs_lat_deg,
            obs_height_m,
            temp_c,
            pressure_mb,
        );
        out.push(ReducedOrbitEntry {
            index: idx + 1,
            epoch: elem.epoch,
            inclination_deg: elem.i,
            asc_node_deg: elem.W,
            arg_perihelion_deg: elem.w,
            semi_major_axis_au: elem.a,
            daily_motion_deg: elem.dm,
            eccentricity: elem.ecc,
            mean_anomaly_deg: elem.M,
            equinox: elem.equinox,
            abs_mag: elem.mag,
            ephemeris: eph,
        });
    }
    out
}

