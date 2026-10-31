//! Keplerian & DE404 ephemeris solver (`kepler.c`, `vearth.c`, and built-in orbits from `aa.c`).

use crate::de404_tables::*;
use crate::earth_orientation::{epsiln, precess};
use crate::gplan::{g3plan, gmoon, gplan, mean_elements};
use crate::math_time::{mod360, modtp, zatan2, DTR, EMRAT, J2000, RTD, TPI};

/// Identifies which DE404 table (if any) a built-in planet uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum De404Body {
    None,
    Mercury,
    Venus,
    EarthBary,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune,
    Pluto,
}

/// Orbital elements structure (`struct orbit` in `kep.h`).
#[derive(Debug, Clone)]
pub struct OrbitElement {
    pub obname: String,
    pub epoch: f64,
    pub i: f64,
    pub W: f64,
    pub w: f64,
    pub a: f64,
    pub dm: f64,
    pub ecc: f64,
    pub M: f64,
    pub equinox: f64,
    pub mag: f64,
    pub sdiam: f64,
    pub L: f64,
    pub r: f64,
    pub plat: f64,
    pub de404: De404Body,
}

/// Built-in planetary orbits from `aa.c` (`mercury` through `pluto`), plus `2060 Chiron`.
pub fn builtin_planets() -> Vec<OrbitElement> {
    vec![
        OrbitElement {
            obname: "Mercury".to_string(),
            epoch: 2446800.5,
            i: 7.0048,
            W: 48.177,
            w: 29.074,
            a: 0.387098,
            dm: 4.092339,
            ecc: 0.205628,
            M: 198.7199,
            equinox: 2446800.5,
            mag: -0.42,
            sdiam: 3.36,
            L: 0.0,
            r: 0.0,
            plat: 0.0,
            de404: De404Body::Mercury,
        },
        OrbitElement {
            obname: "Venus".to_string(),
            epoch: 2446800.5,
            i: 3.3946,
            W: 76.561,
            w: 54.889,
            a: 0.723329,
            dm: 1.602130,
            ecc: 0.006757,
            M: 9.0369,
            equinox: 2446800.5,
            mag: -4.40,
            sdiam: 8.34,
            L: 0.0,
            r: 0.0,
            plat: 0.0,
            de404: De404Body::Venus,
        },
        OrbitElement {
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
        },
        OrbitElement {
            obname: "Mars".to_string(),
            epoch: 2446800.5,
            i: 1.8498,
            W: 49.457,
            w: 286.343,
            a: 1.523712,
            dm: 0.524023,
            ecc: 0.093423,
            M: 23.4699,
            equinox: 2446800.5,
            mag: -1.52,
            sdiam: 4.68,
            L: 0.0,
            r: 0.0,
            plat: 0.0,
            de404: De404Body::Mars,
        },
        OrbitElement {
            obname: "Jupiter".to_string(),
            epoch: 2446800.5,
            i: 1.3051,
            W: 100.358,
            w: 275.129,
            a: 5.20265,
            dm: 0.0830948,
            ecc: 0.048100,
            M: 344.5086,
            equinox: 2446800.5,
            mag: -9.40,
            sdiam: 98.44,
            L: 0.0,
            r: 0.0,
            plat: 0.0,
            de404: De404Body::Jupiter,
        },
        OrbitElement {
            obname: "Saturn".to_string(),
            epoch: 2446800.5,
            i: 2.4858,
            W: 113.555,
            w: 337.969,
            a: 9.54050,
            dm: 0.0334979,
            ecc: 0.052786,
            M: 173.1654,
            equinox: 2446800.5,
            mag: -8.88,
            sdiam: 82.73,
            L: 0.0,
            r: 0.0,
            plat: 0.0,
            de404: De404Body::Saturn,
        },
        OrbitElement {
            obname: "Uranus".to_string(),
            epoch: 2446800.5,
            i: 0.7738,
            W: 73.994,
            w: 98.746,
            a: 19.2233,
            dm: 0.0116943,
            ecc: 0.045682,
            M: 84.8516,
            equinox: 2446800.5,
            mag: -7.19,
            sdiam: 35.02,
            L: 0.0,
            r: 0.0,
            plat: 0.0,
            de404: De404Body::Uranus,
        },
        OrbitElement {
            obname: "Neptune".to_string(),
            epoch: 2446800.5,
            i: 1.7697,
            W: 131.677,
            w: 250.623,
            a: 30.1631,
            dm: 0.0059610,
            ecc: 0.009019,
            M: 254.2568,
            equinox: 2446800.5,
            mag: -6.87,
            sdiam: 33.50,
            L: 0.0,
            r: 0.0,
            plat: 0.0,
            de404: De404Body::Neptune,
        },
        OrbitElement {
            obname: "Pluto".to_string(),
            epoch: 2446640.5,
            i: 17.1346,
            W: 110.204,
            w: 114.21,
            a: 39.4633,
            dm: 0.0039792,
            ecc: 0.248662,
            M: 355.0554,
            equinox: 2446640.5,
            mag: -1.0,
            sdiam: 2.07,
            L: 0.0,
            r: 0.0,
            plat: 0.0,
            de404: De404Body::Pluto,
        },
        OrbitElement {
            obname: "2060 Chiron".to_string(),
            epoch: 2450100.5,
            i: 6.93122,
            W: 209.39100,
            w: 339.55677,
            a: 13.702248,
            dm: 0.0,
            ecc: 0.3822777,
            M: 359.75470,
            equinox: J2000,
            mag: 6.5,
            sdiam: 0.1,
            L: 0.0,
            r: 0.0,
            plat: 0.0,
            de404: De404Body::None,
        },
    ]
}

/// Adjust position from Earth-Moon barycenter to Earth (`embofs` in `kepler.c`).
/// `jed`: Julian date (TDT), `ea`: J2000 equatorial rectangular coordinates of EMB.
/// Returns `(earth_rect_j2000, helio_radius_au)`.
pub fn embofs(jed: f64, ea: &mut [f64; 3]) -> f64 {
    let (mut pm, _polm) = gmoon(jed);
    precess(&mut pm, jed, 1); // Precess geocentric Moon vector from date to J2000

    let mut b = 0.0;
    for i in 0..3 {
        ea[i] -= pm[i] / (EMRAT + 1.0);
        b += ea[i] * ea[i];
    }
    b.sqrt()
}

/// Solve Kepler's equation / DE404 ephemeris (`kepler` in `kepler.c`).
/// Computes heliocentric equatorial J2000 rectangular coordinates `rect` and ecliptic polar coordinates `polar` (`[lon, lat, r]`).
pub fn kepler(jed: f64, e: &mut OrbitElement) -> ([f64; 3], [f64; 3]) {
    let mut rect = [0.0f64; 3];
    let mut polar = [0.0f64; 3];

    let coseps_j2000 = 0.917482062069181825f64;
    let sineps_j2000 = 0.3977771559319137015f64;

    if e.de404 != De404Body::None {
        if e.de404 == De404Body::EarthBary {
            rect = g3plan(jed, &EAR404, 3);
        } else {
            let tbl = match e.de404 {
                De404Body::Mercury => &MER404,
                De404Body::Venus => &VEN404,
                De404Body::Mars => &MAR404,
                De404Body::Jupiter => &JUP404,
                De404Body::Saturn => &SAT404,
                De404Body::Uranus => &URA404,
                De404Body::Neptune => &NEP404,
                De404Body::Pluto => &PLU404,
                _ => unreachable!(),
            };
            rect = gplan(jed, tbl);
        }

        // Compute J2000 ecliptic polar coordinates from J2000 equatorial rectangular
        let x = rect[0];
        let y_ec = coseps_j2000 * rect[1] + sineps_j2000 * rect[2];
        let z_ec = -sineps_j2000 * rect[1] + coseps_j2000 * rect[2];
        let r = (x * x + y_ec * y_ec + z_ec * z_ec).sqrt();
        let lon = zatan2(x, y_ec);
        let lat = (z_ec / r).asin();
        polar = [lon, lat, r];

        e.r = r;
        e.L = lon;
        e.plat = lat;

        if e.de404 == De404Body::EarthBary {
            let r_earth = embofs(jed, &mut rect);
            e.r = r_earth;
            let x = rect[0];
            let y_ec = coseps_j2000 * rect[1] + sineps_j2000 * rect[2];
            let z_ec = -sineps_j2000 * rect[1] + coseps_j2000 * rect[2];
            let lon_e = zatan2(x, y_ec);
            let lat_e = (z_ec / r_earth).asin();
            polar = [lon_e, lat_e, r_earth];
            e.L = lon_e;
            e.plat = lat_e;
        }

        return (rect, polar);
    }

    // Keplerian orbit calculation (`kepler.c`) for asteroids, comets, Chiron, etc.
    let mut alat: f64;
    let mut E = 0.0f64;
    let mut M: f64;
    let mut W: f64;
    let mut temp: f64;

    // Call program to compute position of Earth if this were Earth, else Keplerian
    let epoch = e.epoch;
    let mut inclination = e.i;
    let mut ascnode = e.W * DTR;
    let argperih = e.w;
    let meandistance = e.a;
    let mut dailymotion = e.dm;
    let eccent = e.ecc;
    let meananomaly = e.M;

    if dailymotion == 0.0 {
        if eccent > 1.0 {
            dailymotion = 0.9856076686 / (e.a * (e.a * (eccent - 1.0)).sqrt());
        } else {
            dailymotion = 0.9856076686 / (e.a * e.a.sqrt());
        }
    }

    if eccent == 1.0 {
        // Parabolic orbit (`epoch` is perihelion date, `a` is q)
        let k_gauss = 0.01720209895;
        let q = meandistance;
        let d_time = jed - epoch;
        let w_para = (3.0 * k_gauss / (2.0f64).sqrt()) * d_time / (q * q.sqrt());
        let y_para = (w_para * 0.5 + (w_para * w_para * 0.25 + 1.0).sqrt()).cbrt();
        let s_true = y_para - 1.0 / y_para;
        let v_true = 2.0 * s_true.atan();
        let r = q * (1.0 + s_true * s_true);
        alat = v_true + DTR * argperih;
        e.r = r;
    } else if eccent > 1.0 {
        // Hyperbolic orbit
        let meandistance_h = meandistance / (eccent - 1.0);
        let dailymotion_h = 0.9856076686 / (meandistance_h * meandistance_h.sqrt());
        M = dailymotion_h * (jed - epoch);
        E = M * DTR;
        for _ in 0..25 {
            temp = E - (eccent * E.sinh() - E - M * DTR) / (eccent * E.cosh() - 1.0);
            if (temp - E).abs() <= 1.0e-11 {
                E = temp;
                break;
            }
            E = temp;
        }
        let r = meandistance_h * (eccent * E.cosh() - 1.0);
        e.r = r;
        let num = (eccent + 1.0).sqrt() * (0.5 * E).sinh();
        let den = (eccent - 1.0).sqrt() * (0.5 * E).cosh();
        let v_true = 2.0 * num.atan2(den);
        alat = v_true + DTR * argperih;
    } else {
        // Elliptical orbit
        M = dailymotion * (jed - epoch) + meananomaly;
        M = mod360(M);
        M *= DTR;

        E = M + eccent * M.sin() + 0.5 * eccent * eccent * (2.0 * M).sin();
        for _ in 0..25 {
            temp = E - (E - eccent * E.sin() - M) / (1.0 - eccent * E.cos());
            if (temp - E).abs() <= 1.0e-11 {
                E = temp;
                break;
            }
            E = temp;
        }
        let r = meandistance * (1.0 - eccent * E.cos());
        e.r = r;
        let num = (1.0 + eccent).sqrt() * (0.5 * E).sin();
        let den = (1.0 - eccent).sqrt() * (0.5 * E).cos();
        let v_true = 2.0 * zatan2(den, num);
        alat = v_true + DTR * argperih;
    }

    // From heliocentric orbital plane to ecliptic coordinates of `e.equinox`
    inclination *= DTR;
    let salp = alat.sin();
    let calp = alat.cos();
    let sbeta = salp * inclination.sin();
    let cbeta = (1.0 - sbeta * sbeta).sqrt();

    let mut l_orb = zatan2(calp, salp * inclination.cos()) + ascnode;
    l_orb = modtp(l_orb);
    let b_orb = sbeta.asin();

    polar[0] = l_orb;
    polar[1] = b_orb;
    polar[2] = e.r;

    e.L = polar[0];
    e.plat = polar[1];

    // Convert ecliptic polar (at `e.equinox`) to equatorial rectangular (at `e.equinox`)
    let eps_eq = if (e.equinox - J2000).abs() < 1.0e-4 {
        epsiln(J2000)
    } else {
        epsiln(e.equinox)
    };
    let coseps = eps_eq.cos();
    let sineps = eps_eq.sin();

    let x_ec = e.r * cbeta * l_orb.cos();
    let y_ec = e.r * cbeta * l_orb.sin();
    let z_ec = e.r * sbeta;

    rect[0] = x_ec;
    rect[1] = coseps * y_ec - sineps * z_ec;
    rect[2] = sineps * y_ec + coseps * z_ec;

    // Precess from `e.equinox` to J2000 if needed
    if (e.equinox - J2000).abs() > 1.0e-4 {
        precess(&mut rect, e.equinox, 1);
    }

    (rect, polar)
}

/// Calculate velocity vector of the Earth (`vearth.c`) in AU/day (J2000 equatorial rectangular).
/// Uses finite difference over `t = 0.005` days exactly as in `vearth.c`.
pub fn velearth(jed: f64) -> [f64; 3] {
    let dt = 0.005;
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
    let (p_minus, _) = kepler(jed - dt, &mut earth_el);
    let (p_plus, _) = kepler(jed, &mut earth_el);
    [
        (p_plus[0] - p_minus[0]) / dt,
        (p_plus[1] - p_minus[1]) / dt,
        (p_plus[2] - p_minus[2]) / dt,
    ]
}

/// Osculating Lunar Perigee (True Black Moon Lilith) and True Lunar Ascending Node (`mean_elements` + osculating series).
/// Returns `(mean_lilith_rad, true_lilith_rad, mean_node_rad, true_node_rad)` re ecliptic of date.
pub fn compute_lunar_nodes_and_lilith(jed: f64) -> (f64, f64, f64, f64) {
    let me = mean_elements(jed);
    let d = me.mean_elongation_rad;
    let m_sun = me.args[11];
    let m_moon = me.moon_mean_anomaly_rad;
    let f = me.args[10];

    // True (osculating) lunar ascending node perturbation (Meeus / Chapront)
    let d_node_deg = -1.4979 * (2.0 * (d - f)).sin()
        - 0.1500 * m_sun.sin()
        - 0.1226 * (2.0 * d).sin()
        + 0.1176 * (2.0 * f).sin()
        - 0.0801 * (2.0 * (m_moon - f)).sin();
    let true_node_rad = modtp(me.lunar_node_rad + d_node_deg * DTR);

    // True (osculating) lunar perigee / True Black Moon Lilith perturbation
    let d_perigee_deg = 1.794 * (2.0 * d - m_moon).sin()
        - 0.406 * (2.0 * d).sin()
        - 0.276 * (2.0 * d - 2.0 * m_moon).sin()
        + 0.196 * (4.0 * d - m_moon).sin()
        - 0.193 * m_moon.sin();
    let true_lilith_rad = modtp(me.lunar_perigee_rad + d_perigee_deg * DTR);

    (
        modtp(me.lunar_perigee_rad),
        true_lilith_rad,
        modtp(me.lunar_node_rad),
        true_node_rad,
    )
}

