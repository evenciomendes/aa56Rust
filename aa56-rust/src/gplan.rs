//! Trigonometric series evaluator for DE404 planetary and lunar tables (`gplan.c`).
//! Implements `mean_elements`, `gplan`, `g3plan`, `g2plan`, `g1plan`, and `gmoon`.

use crate::de404_tables::{PlanTblF64, PlanTblI32, MOONLAT, MOONLR};
use crate::earth_orientation::epsiln;
use crate::math_time::{mods3600, DTR, J2000, RTS, STR, TPI};

const COSEPS_J2000: f64 = 0.917482062069181825;
const SINEPS_J2000: f64 = 0.3977771559319137015;

/// Polynomial coefficients for fundamental mean elements (`gplan.c`).
static SIMON_FREQS: [f64; 9] = [
    538101628.688982,
    210664136.433548,
    129597742.283429,
    68905077.493988,
    10925660.377991,
    4399609.855372,
    1542481.193933,
    786550.320744,
    522722.451795,
];

static SIMON_PHASES: [f64; 9] = [
    252.25090552 * 3600.0,
    181.97980085 * 3600.0,
    100.46645683 * 3600.0,
    355.43299958 * 3600.0,
    34.35151874 * 3600.0,
    50.07744430 * 3600.0,
    314.05500511 * 3600.0,
    304.34866548 * 3600.0,
    860492.1546,
];

#[derive(Debug, Clone, Copy)]
pub struct MeanElements {
    pub args: [f64; 18],
    pub ea_arcsec: f64,
    pub lp_equinox: f64,
    pub nf_arcsec: f64,
    /// Mean longitude of Lunar Perigee (Black Moon Lilith) in radians (of date)
    pub lunar_perigee_rad: f64,
    /// Mean longitude of Lunar Ascending Node (Rahu) in radians (of date)
    pub lunar_node_rad: f64,
    /// Mean elongation D in radians
    pub mean_elongation_rad: f64,
    /// Mean anomaly of Moon l in radians
    pub moon_mean_anomaly_rad: f64,
}

/// Evaluate fundamental mean elements of the solar system and Moon at Julian date `jed` (`mean_elements` in `gplan.c`).
pub fn mean_elements(jed: f64) -> MeanElements {
    let mut args = [0.0f64; 18];
    let t = (jed - J2000) / 36525.0;
    let t2 = t * t;

    // Mean longitudes of planets (Simon et al., 1994)
    for i in 0..9 {
        let x = mods3600(SIMON_FREQS[i] * t + SIMON_PHASES[i]);
        args[i] = x * STR;
    }

    // 9: Mean elongation of Moon from Sun (D)
    let mut x = 1.6029616019939659e+09 * t + 1.0722612202445078e+06;
    x += (((((-3.207663637426e-013 * t + 2.555243317839e-011) * t + 2.560078201452e-009)
        * t
        - 3.702060118571e-005)
        * t
        + 6.9492746836058421e-03)
        * t
        - 6.7352202374457519e+00)
        * t2;
    let d_arcsec = mods3600(x);
    args[9] = STR * d_arcsec;

    // 10: Mean distance of Moon from its ascending node (F)
    x = 1.7395272628437717e+09 * t + 3.3577951412884740e+05;
    x += (((((4.474984866301e-013 * t + 4.189032191814e-011) * t - 2.790392351314e-009) * t
        - 2.165750777942e-006)
        * t
        - 7.5311878482337989e-04)
        * t
        - 1.3117809789650071e+01)
        * t2;
    let nf_arcsec = mods3600(x);
    args[10] = STR * nf_arcsec;

    // 11: Mean anomaly of Sun (l')
    x = 1.2959658102304320e+08 * t + 1.2871027407441526e+06;
    x += ((((((((-1.62e-20 * t - 1.0390e-17) * t - 3.83508e-15) * t + 4.237343e-13) * t
        + 8.8555011e-11)
        * t
        - 4.77258489e-8)
        * t
        - 1.1297037031e-5)
        * t
        + 8.7473717367324703e-05)
        * t
        - 5.5281306421783094e-01)
        * t2;
    args[11] = STR * mods3600(x);

    // 12: Mean anomaly of Moon (l)
    x = 1.7179159228846793e+09 * t + 4.8586817465825332e+05;
    x += (((((-1.755312760154e-012 * t + 3.452144225877e-011) * t - 2.506365935364e-008)
        * t
        - 2.536291235258e-004)
        * t
        + 5.2099641302735818e-02)
        * t
        + 3.1501359071894147e+01)
        * t2;
    let l_moon_arcsec = mods3600(x);
    args[12] = STR * l_moon_arcsec;

    // 13: Mean longitude of Moon, re mean ecliptic and equinox of date (L)
    x = 1.7325643720442266e+09 * t + 7.8593980921052420e+05;
    x += (((((-7.200592540556e-014 * t + 2.235210987108e-010) * t - 1.024222633731e-008)
        * t
        - 6.073960534117e-005)
        * t
        + 6.9017248528380490e-03)
        * t
        - 5.6550460027471399e+00)
        * t2;
    let lp_equinox = mods3600(x);
    args[13] = STR * lp_equinox;

    // 14: Precession of the equinox
    x = (((((((((-8.66e-20 * t - 4.759e-17) * t + 2.424e-15) * t + 1.3095e-12) * t
        + 1.7451e-10)
        * t
        - 1.8055e-8)
        * t
        - 0.0000235316)
        * t
        + 0.000076)
        * t
        + 1.105414)
        * t
        + 5028.791959)
        * t;
    // Moon's longitude re fixed J2000 equinox
    let ea_arcsec = args[2] * RTS;
    args[14] = (lp_equinox - x) * STR;

    // 15: Free librations / lunar node
    // Mean longitude of lunar ascending node Omega = L - F
    let node_arcsec = mods3600(lp_equinox - nf_arcsec);
    // Mean longitude of lunar perigee Lp = L - l
    let perigee_arcsec = mods3600(lp_equinox - l_moon_arcsec);

    x = (4.0729050039324226e+07 * t + 1.2626554492184688e+06) - x;
    args[15] = STR * mods3600(x);
    args[16] = STR * mods3600(5.36486787828342e+06 * t + 9.3308862598336e+05);
    args[17] = STR * node_arcsec;

    MeanElements {
        args,
        ea_arcsec,
        lp_equinox,
        nf_arcsec,
        lunar_perigee_rad: STR * perigee_arcsec,
        lunar_node_rad: STR * node_arcsec,
        mean_elongation_rad: STR * d_arcsec,
        moon_mean_anomaly_rad: STR * l_moon_arcsec,
    }
}

/// Generic program to accumulate sum of trigonometric series for 3 variables (`gplan` in `gplan.c`).
/// Used for `MER404`, `VEN404`, `MAR404`, `JUP404`, `SAT404`, `URA404`, `NEP404`, `PLU404` (`f64` tables).
/// Outputs J2000 heliocentric equatorial rectangular coordinates in `pobj`.
pub fn gplan(jed: f64, plan: &PlanTblF64) -> [f64; 3] {
    let me = mean_elements(jed);
    let t = (jed - J2000) / plan.timescale;
    let n = plan.maxargs as usize;

    let mut ss = [[0.0f64; 41]; 18];
    let mut cc = [[0.0f64; 41]; 18];

    for i in 0..n {
        let j = plan.max_harmonic[i] as usize;
        if j > 0 {
            let sr = me.args[i].sin();
            let cr = me.args[i].cos();
            ss[i][0] = sr;
            cc[i][0] = cr;
            for k in 1..j {
                ss[i][k] = cr * ss[i][k - 1] + sr * cc[i][k - 1];
                cc[i][k] = cr * cc[i][k - 1] - sr * ss[i][k - 1];
            }
        }
    }

    let mut sl = 0.0f64;
    let mut sb = 0.0f64;
    let mut sr = 0.0f64;

    let mut p_idx = 0usize;
    let mut pl_idx = 0usize;
    let mut pb_idx = 0usize;
    let mut pr_idx = 0usize;

    loop {
        let m = plan.arg_tbl[p_idx];
        p_idx += 1;
        if m < 0 {
            break;
        }
        let mut su = 0.0f64;
        let mut cu = 1.0f64;

        for _ in 0..m {
            let mut j = plan.arg_tbl[p_idx] as i32;
            p_idx += 1;
            let k1 = (plan.arg_tbl[p_idx] as usize) - 1;
            p_idx += 1;

            if j != 0 {
                let mut k = j.unsigned_abs() as usize - 1;
                let mut s = ss[k1][k];
                if j < 0 {
                    s = -s;
                }
                let c = cc[k1][k];
                let k2 = cu * s + su * c;
                cu = cu * c - su * s;
                su = k2;
            }
        }

        let nt = plan.arg_tbl[p_idx] as i32;
        p_idx += 1;

        if nt == 0 {
            let cv = plan.lon_tbl[pl_idx];
            let sv = plan.lon_tbl[pl_idx + 1];
            pl_idx += 2;
            sl += cu * cv + su * sv;

            let cv = plan.lat_tbl[pb_idx];
            let sv = plan.lat_tbl[pb_idx + 1];
            pb_idx += 2;
            sb += cu * cv + su * sv;

            let cv = plan.rad_tbl[pr_idx];
            let sv = plan.rad_tbl[pr_idx + 1];
            pr_idx += 2;
            sr += cu * cv + su * sv;
        } else {
            let mut k1 = 0.0f64;
            let mut k2 = 0.0f64;
            let mut kb = 0.0f64;
            let mut kc = 0.0f64;
            let mut kr = 0.0f64;
            let mut ks = 0.0f64;
            for _ in 0..nt {
                let cv = plan.lon_tbl[pl_idx];
                let sv = plan.lon_tbl[pl_idx + 1];
                pl_idx += 2;
                k1 = k1 * t + cv;
                k2 = k2 * t + sv;

                let cv = plan.lat_tbl[pb_idx];
                let sv = plan.lat_tbl[pb_idx + 1];
                pb_idx += 2;
                kb = kb * t + cv;
                kc = kc * t + sv;

                let cv = plan.rad_tbl[pr_idx];
                let sv = plan.rad_tbl[pr_idx + 1];
                pr_idx += 2;
                kr = kr * t + cv;
                ks = ks * t + sv;
            }
            if m != 0 {
                sl += cu * k1 + su * k2;
                sb += cu * kb + su * kc;
                sr += cu * kr + su * ks;
            } else {
                sl += k1;
                sb += kb;
                sr += kr;
            }
        }
    }

    let lon = STR * sl;
    let lat = STR * sb;
    let r = STR * sr * plan.distance + plan.distance;

    // Convert heliocentric ecliptic polar J2000 to heliocentric equatorial rectangular J2000
    let cl = lon.cos();
    let sl_sin = lon.sin();
    let cb = lat.cos();
    let sb_sin = lat.sin();

    let x = r * cb * cl;
    let y = r * cb * sl_sin;
    let z = r * sb_sin;

    [
        x,
        COSEPS_J2000 * y - SINEPS_J2000 * z,
        SINEPS_J2000 * y + COSEPS_J2000 * z,
    ]
}

/// Accumulate sum of trigonometric series for 3 variables (`g3plan` in `gplan.c`).
/// Used for `EAR404` (`i32` tables, `objnum = 3`).
/// Outputs J2000 heliocentric equatorial rectangular coordinates of the Earth-Moon barycenter.
pub fn g3plan(jed: f64, plan: &PlanTblI32, objnum: usize) -> [f64; 3] {
    let me = mean_elements(jed);
    let t = (jed - J2000) / plan.timescale;
    let n = plan.maxargs as usize;

    let mut ss = [[0.0f64; 41]; 18];
    let mut cc = [[0.0f64; 41]; 18];

    for i in 0..n {
        let j = plan.max_harmonic[i] as usize;
        if j > 0 {
            let u = if objnum == 3 && i == 2 {
                me.ea_arcsec * STR
            } else {
                me.args[i]
            };
            let sr = u.sin();
            let cr = u.cos();
            ss[i][0] = sr;
            cc[i][0] = cr;
            for k in 1..j {
                ss[i][k] = cr * ss[i][k - 1] + sr * cc[i][k - 1];
                cc[i][k] = cr * cc[i][k - 1] - sr * ss[i][k - 1];
            }
        }
    }

    let mut sl = 0.0f64;
    let mut sb = 0.0f64;
    let mut sr = 0.0f64;

    let mut p_idx = 0usize;
    let mut pl_idx = 0usize;
    let mut pb_idx = 0usize;
    let mut pr_idx = 0usize;

    loop {
        let m = plan.arg_tbl[p_idx];
        p_idx += 1;
        if m < 0 {
            break;
        }
        let mut su = 0.0f64;
        let mut cu = 1.0f64;

        for _ in 0..m {
            let j = plan.arg_tbl[p_idx] as i32;
            p_idx += 1;
            let k1 = (plan.arg_tbl[p_idx] as usize) - 1;
            p_idx += 1;

            if j != 0 {
                let k = j.unsigned_abs() as usize - 1;
                let mut s = ss[k1][k];
                if j < 0 {
                    s = -s;
                }
                let c = cc[k1][k];
                let k2 = cu * s + su * c;
                cu = cu * c - su * s;
                su = k2;
            }
        }

        let mut k = plan.arg_tbl[p_idx] as i32;
        p_idx += 1;

        let mut su1 = plan.lon_tbl[pl_idx] as f64;
        let mut cu1 = plan.lon_tbl[pl_idx + 1] as f64;
        pl_idx += 2;

        let mut sb1 = plan.lat_tbl[pb_idx] as f64;
        let mut cb1 = plan.lat_tbl[pb_idx + 1] as f64;
        pb_idx += 2;

        let mut sr1 = plan.rad_tbl[pr_idx] as f64;
        let mut cr1 = plan.rad_tbl[pr_idx + 1] as f64;
        pr_idx += 2;

        if k > 0 {
            let mut lt = t;
            while k > 0 {
                su1 += lt * (plan.lon_tbl[pl_idx] as f64);
                cu1 += lt * (plan.lon_tbl[pl_idx + 1] as f64);
                pl_idx += 2;

                sb1 += lt * (plan.lat_tbl[pb_idx] as f64);
                cb1 += lt * (plan.lat_tbl[pb_idx + 1] as f64);
                pb_idx += 2;

                sr1 += lt * (plan.rad_tbl[pr_idx] as f64);
                cr1 += lt * (plan.rad_tbl[pr_idx + 1] as f64);
                pr_idx += 2;

                lt *= t;
                k -= 1;
            }
        }

        sl += cu1 * cu + su1 * su;
        sb += cb1 * cu + sb1 * su;
        sr += cr1 * cu + sr1 * su;
    }

    let mut lon_sec = plan.trunclvl * sl;
    if objnum == 3 {
        lon_sec += me.ea_arcsec;
    }
    let lon = STR * lon_sec;
    let lat = STR * plan.trunclvl * sb;
    let r = plan.distance * (1.0 + STR * plan.trunclvl * sr);

    let cl = lon.cos();
    let sl_sin = lon.sin();
    let cb = lat.cos();
    let sb_sin = lat.sin();

    let x = r * cb * cl;
    let y = r * cb * sl_sin;
    let z = r * sb_sin;

    [
        x,
        COSEPS_J2000 * y - SINEPS_J2000 * z,
        SINEPS_J2000 * y + COSEPS_J2000 * z,
    ]
}

/// Accumulate sum of trigonometric series for 2 variables (`g2plan` in `gplan.c`).
/// Used for `MOONLR` (longitude and radius perturbations in arcseconds).
pub fn g2plan(jed: f64, plan: &PlanTblI32, me: &MeanElements) -> (f64, f64) {
    let t = (jed - J2000) / plan.timescale;
    let n = plan.maxargs as usize;

    let mut ss = [[0.0f64; 41]; 18];
    let mut cc = [[0.0f64; 41]; 18];

    for i in 0..n {
        let j = plan.max_harmonic[i] as usize;
        if j > 0 {
            let sr = me.args[i].sin();
            let cr = me.args[i].cos();
            ss[i][0] = sr;
            cc[i][0] = cr;
            for k in 1..j {
                ss[i][k] = cr * ss[i][k - 1] + sr * cc[i][k - 1];
                cc[i][k] = cr * cc[i][k - 1] - sr * ss[i][k - 1];
            }
        }
    }

    let mut sl = 0.0f64;
    let mut sr = 0.0f64;

    let mut p_idx = 0usize;
    let mut pl_idx = 0usize;
    let mut pr_idx = 0usize;

    loop {
        let m = plan.arg_tbl[p_idx];
        p_idx += 1;
        if m < 0 {
            break;
        }
        let mut su = 0.0f64;
        let mut cu = 1.0f64;

        for _ in 0..m {
            let j = plan.arg_tbl[p_idx] as i32;
            p_idx += 1;
            let k1 = (plan.arg_tbl[p_idx] as usize) - 1;
            p_idx += 1;

            if j != 0 {
                let k = j.unsigned_abs() as usize - 1;
                let mut s = ss[k1][k];
                if j < 0 {
                    s = -s;
                }
                let c = cc[k1][k];
                let k2 = cu * s + su * c;
                cu = cu * c - su * s;
                su = k2;
            }
        }

        let mut k = plan.arg_tbl[p_idx] as i32;
        p_idx += 1;

        let mut su1 = plan.lon_tbl[pl_idx] as f64;
        let mut cu1 = plan.lon_tbl[pl_idx + 1] as f64;
        pl_idx += 2;

        let mut sr1 = plan.rad_tbl[pr_idx] as f64;
        let mut cr1 = plan.rad_tbl[pr_idx + 1] as f64;
        pr_idx += 2;

        if k > 0 {
            let mut lt = t;
            while k > 0 {
                su1 += lt * (plan.lon_tbl[pl_idx] as f64);
                cu1 += lt * (plan.lon_tbl[pl_idx + 1] as f64);
                pl_idx += 2;

                sr1 += lt * (plan.rad_tbl[pr_idx] as f64);
                cr1 += lt * (plan.rad_tbl[pr_idx + 1] as f64);
                pr_idx += 2;

                lt *= t;
                k -= 1;
            }
        }

        sl += cu1 * cu + su1 * su;
        sr += cr1 * cu + sr1 * su;
    }

    (plan.trunclvl * sl, plan.trunclvl * sr)
}

/// Accumulate sum of trigonometric series for 1 variable (`g1plan` in `gplan.c`).
/// Used for `MOONLAT` (latitude in arcseconds).
pub fn g1plan(jed: f64, plan: &PlanTblI32, me: &MeanElements) -> f64 {
    let t = (jed - J2000) / plan.timescale;
    let n = plan.maxargs as usize;

    let mut ss = [[0.0f64; 41]; 18];
    let mut cc = [[0.0f64; 41]; 18];

    for i in 0..n {
        let j = plan.max_harmonic[i] as usize;
        if j > 0 {
            let sr = me.args[i].sin();
            let cr = me.args[i].cos();
            ss[i][0] = sr;
            cc[i][0] = cr;
            for k in 1..j {
                ss[i][k] = cr * ss[i][k - 1] + sr * cc[i][k - 1];
                cc[i][k] = cr * cc[i][k - 1] - sr * ss[i][k - 1];
            }
        }
    }

    let mut sl = 0.0f64;
    let mut p_idx = 0usize;
    let mut pl_idx = 0usize;

    loop {
        let m = plan.arg_tbl[p_idx];
        p_idx += 1;
        if m < 0 {
            break;
        }
        let mut su = 0.0f64;
        let mut cu = 1.0f64;

        for _ in 0..m {
            let j = plan.arg_tbl[p_idx] as i32;
            p_idx += 1;
            let k1 = (plan.arg_tbl[p_idx] as usize) - 1;
            p_idx += 1;

            if j != 0 {
                let k = j.unsigned_abs() as usize - 1;
                let mut s = ss[k1][k];
                if j < 0 {
                    s = -s;
                }
                let c = cc[k1][k];
                let k2 = cu * s + su * c;
                cu = cu * c - su * s;
                su = k2;
            }
        }

        let mut k = plan.arg_tbl[p_idx] as i32;
        p_idx += 1;

        let mut su1 = plan.lon_tbl[pl_idx] as f64;
        let mut cu1 = plan.lon_tbl[pl_idx + 1] as f64;
        pl_idx += 2;

        if k > 0 {
            let mut lt = t;
            while k > 0 {
                su1 += lt * (plan.lon_tbl[pl_idx] as f64);
                cu1 += lt * (plan.lon_tbl[pl_idx + 1] as f64);
                pl_idx += 2;
                lt *= t;
                k -= 1;
            }
        }

        sl += cu1 * cu + su1 * su;
    }

    plan.trunclvl * sl
}

/// Compute geocentric Moon coordinates (`gmoon` in `gplan.c`).
/// Returns `(rect_of_date, pol_of_date)` where `pol` is `[lon_rad, lat_rad, dist_au]` re mean ecliptic and equinox of date,
/// and `rect` is geocentric equatorial rectangular coordinates of date in AU.
pub fn gmoon(jed: f64) -> ([f64; 3], [f64; 3]) {
    let me = mean_elements(jed);
    let (dlon_arcsec, drad_arcsec) = g2plan(jed, &MOONLR, &me);
    let mut x = dlon_arcsec + me.lp_equinox;
    if x < -6.48e5 {
        x += 1.296e6;
    }
    if x > 6.48e5 {
        x -= 1.296e6;
    }
    let lon = STR * x;
    let lat_arcsec = g1plan(jed, &MOONLAT, &me);
    let lat = STR * lat_arcsec;
    let r = (1.0 + STR * drad_arcsec) * MOONLR.distance;
    let pol = [lon, lat, r];

    // Convert ecliptic polar of date to equatorial rectangular of date (`epsiln(J)`)
    let eps = epsiln(jed);
    let coseps = eps.cos();
    let sineps = eps.sin();

    let cosB = lat.cos();
    let sinB = lat.sin();
    let cosL = lon.cos();
    let sinL = lon.sin();

    let x_ec = cosB * cosL * r;
    let y_ec = cosB * sinL * r;
    let z_ec = sinB * r;

    let rect = [
        x_ec,
        coseps * y_ec - sineps * z_ec,
        sineps * y_ec + coseps * z_ec,
    ];
    (rect, pol)
}

