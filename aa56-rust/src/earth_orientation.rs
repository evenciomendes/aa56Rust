//! Precessão (`precess.c` WILLIAMS/DE403), Nutação IAU 1980 completa (`nutate.c`),
//! Obliquidade da Eclíptica (`epsiln.c` WILLIAMS/DE403), Tempo Sideral Aparente (`sidrlt.c` DE403),
//! Refração (`refrac.c`) e Redução Topocêntrica (`altaz.c`, `diurab.c`, `diurpx.c`, `kfiles.c`).

use crate::math_time::{mod3600, zatan2, AEARTH_M, AU_KM, DTR, FLAT_INV, J2000, RTD, STR};

/// Obliquidade média da eclíptica na data Juliana `jd` (`epsiln.c` - fórmula WILLIAMS = 1 / DE403).
/// Retorna `(eps, cos_eps, sin_eps)` em radianos.
pub fn epsiln(jd: f64) -> (f64, f64, f64) {
    let t = (jd - J2000) / 365250.0; // T /= 10.0 em epsiln.c
    let eps_arcsec = (((((((((2.45e-10 * t + 5.79e-9) * t + 2.787e-7) * t + 7.12e-7) * t
        - 3.905e-5)
        * t
        - 2.4967e-3)
        * t
        - 5.138e-3)
        * t
        + 1.9989)
        * t
        - 0.0175)
        * t
        - 468.33960)
        * t
        + 84381.406173;
    let eps = eps_arcsec * STR;
    (eps, eps.cos(), eps.sin())
}

/// Tabela completa de 105 linhas (9 colunas) da série de nutação IAU 1980 (`nutate.c`).
const NT: [i16; 105 * 9] = [
     0, 0, 0, 0, 2, 2062, 2,-895, 5,
    -2, 0, 2, 0, 1, 46, 0,-24, 0,
     2, 0,-2, 0, 0, 11, 0, 0, 0,
    -2, 0, 2, 0, 2,-3, 0, 1, 0,
     1,-1, 0,-1, 0,-3, 0, 0, 0,
     0,-2, 2,-2, 1,-2, 0, 1, 0,
     2, 0,-2, 0, 1, 1, 0, 0, 0,
     0, 0, 2,-2, 2,-13187,-16, 5736,-31,
     0, 1, 0, 0, 0, 1426,-34, 54,-1,
     0, 1, 2,-2, 2,-517, 12, 224,-6,
     0,-1, 2,-2, 2, 217,-5,-95, 3,
     0, 0, 2,-2, 1, 129, 1,-70, 0,
     2, 0, 0,-2, 0, 48, 0, 1, 0,
     0, 0, 2,-2, 0,-22, 0, 0, 0,
     0, 2, 0, 0, 0, 17,-1, 0, 0,
     0, 1, 0, 0, 1,-15, 0, 9, 0,
     0, 2, 2,-2, 2,-16, 1, 7, 0,
     0,-1, 0, 0, 1,-12, 0, 6, 0,
    -2, 0, 0, 2, 1,-6, 0, 3, 0,
     0,-1, 2,-2, 1,-5, 0, 3, 0,
     2, 0, 0,-2, 1, 4, 0,-2, 0,
     0, 1, 2,-2, 1, 4, 0,-2, 0,
     1, 0, 0,-1, 0,-4, 0, 0, 0,
     2, 1, 0,-2, 0, 1, 0, 0, 0,
     0, 0,-2, 2, 1, 1, 0, 0, 0,
     0, 1,-2, 2, 0,-1, 0, 0, 0,
     0, 1, 0, 0, 2, 1, 0, 0, 0,
    -1, 0, 0, 1, 1, 1, 0, 0, 0,
     0, 1, 2,-2, 0,-1, 0, 0, 0,
     0, 0, 2, 0, 2,-2274,-2, 977,-5,
     1, 0, 0, 0, 0, 712, 1,-7, 0,
     0, 0, 2, 0, 1,-386,-4, 200, 0,
     1, 0, 2, 0, 2,-301, 0, 129,-1,
     1, 0, 0,-2, 0,-158, 0,-1, 0,
    -1, 0, 2, 0, 2, 123, 0,-53, 0,
     0, 0, 0, 2, 0, 63, 0,-2, 0,
     1, 0, 0, 0, 1, 63, 1,-33, 0,
    -1, 0, 0, 0, 1,-58,-1, 32, 0,
    -1, 0, 2, 2, 2,-59, 0, 26, 0,
     1, 0, 2, 0, 1,-51, 0, 27, 0,
     0, 0, 2, 2, 2,-38, 0, 16, 0,
     2, 0, 0, 0, 0, 29, 0,-1, 0,
     1, 0, 2,-2, 2, 29, 0,-12, 0,
     2, 0, 2, 0, 2,-31, 0, 13, 0,
     0, 0, 2, 0, 0, 26, 0,-1, 0,
    -1, 0, 2, 0, 1, 21, 0,-10, 0,
    -1, 0, 0, 2, 1, 16, 0,-8, 0,
     1, 0, 0,-2, 1,-13, 0, 7, 0,
    -1, 0, 2, 2, 1,-10, 0, 5, 0,
     1, 1, 0,-2, 0,-7, 0, 0, 0,
     0, 1, 2, 0, 2, 7, 0,-3, 0,
     0,-1, 2, 0, 2,-7, 0, 3, 0,
     1, 0, 2, 2, 2,-8, 0, 3, 0,
     1, 0, 0, 2, 0, 6, 0, 0, 0,
     2, 0, 2,-2, 2, 6, 0,-3, 0,
     0, 0, 0, 2, 1,-6, 0, 3, 0,
     0, 0, 2, 2, 1,-7, 0, 3, 0,
     1, 0, 2,-2, 1, 6, 0,-3, 0,
     0, 0, 0,-2, 1,-5, 0, 3, 0,
     1,-1, 0, 0, 0, 5, 0, 0, 0,
     2, 0, 2, 0, 1,-5, 0, 3, 0,
     0, 1, 0,-2, 0,-4, 0, 0, 0,
     1, 0,-2, 0, 0, 4, 0, 0, 0,
     0, 0, 0, 1, 0,-4, 0, 0, 0,
     1, 1, 0, 0, 0,-3, 0, 0, 0,
     1, 0, 2, 0, 0, 3, 0, 0, 0,
     1,-1, 2, 0, 2,-3, 0, 1, 0,
    -1,-1, 2, 2, 2,-3, 0, 1, 0,
    -2, 0, 0, 0, 1,-2, 0, 1, 0,
     3, 0, 2, 0, 2,-3, 0, 1, 0,
     0,-1, 2, 2, 2,-3, 0, 1, 0,
     1, 1, 2, 0, 2, 2, 0,-1, 0,
    -1, 0, 2,-2, 1,-2, 0, 1, 0,
     2, 0, 0, 0, 1, 2, 0,-1, 0,
     1, 0, 0, 0, 2,-2, 0, 1, 0,
     3, 0, 0, 0, 0, 2, 0, 0, 0,
     0, 0, 2, 1, 2, 2, 0,-1, 0,
    -1, 0, 0, 0, 2, 1, 0,-1, 0,
     1, 0, 0,-4, 0,-1, 0, 0, 0,
    -2, 0, 2, 2, 2, 1, 0,-1, 0,
    -1, 0, 2, 4, 2,-2, 0, 1, 0,
     2, 0, 0,-4, 0,-1, 0, 0, 0,
     1, 1, 2,-2, 2, 1, 0,-1, 0,
     1, 0, 2, 2, 1,-1, 0, 1, 0,
    -2, 0, 2, 4, 2,-1, 0, 1, 0,
    -1, 0, 4, 0, 2, 1, 0, 0, 0,
     1,-1, 0,-2, 0, 1, 0, 0, 0,
     2, 0, 2,-2, 1, 1, 0,-1, 0,
     2, 0, 2, 2, 2,-1, 0, 0, 0,
     1, 0, 0, 2, 1,-1, 0, 0, 0,
     0, 0, 4,-2, 2, 1, 0, 0, 0,
     3, 0, 2,-2, 2, 1, 0, 0, 0,
     1, 0, 2,-2, 0,-1, 0, 0, 0,
     0, 1, 2, 0, 1, 1, 0, 0, 0,
    -1,-1, 0, 2, 1, 1, 0, 0, 0,
     0, 0,-2, 0, 1,-1, 0, 0, 0,
     0, 0, 2,-1, 2,-1, 0, 0, 0,
     0, 1, 0, 2, 0,-1, 0, 0, 0,
     1, 0,-2,-2, 0,-1, 0, 0, 0,
     0,-1, 2, 0, 1,-1, 0, 0, 0,
     1, 1, 0,-2, 1,-1, 0, 0, 0,
     1, 0,-2, 2, 0,-1, 0, 0, 0,
     2, 0, 0, 2, 0, 1, 0, 0, 0,
     0, 0, 2, 4, 2,-1, 0, 0, 0,
     0, 1, 0, 1, 0, 1, 0, 0, 0,
];

fn fill_sscc(ss: &mut [[f64; 8]; 5], cc: &mut [[f64; 8]; 5], k: usize, arg: f64, n: usize) {
    let su = arg.sin();
    let cu = arg.cos();
    ss[k][0] = su;
    cc[k][0] = cu;
    let mut sv = 2.0 * su * cu;
    let mut cv = cu * cu - su * su;
    ss[k][1] = sv;
    cc[k][1] = cv;
    for i in 2..n {
        let s = su * cv + cu * sv;
        cv = cu * cv - su * sv;
        sv = s;
        ss[k][i] = sv;
        cc[k][i] = cv;
    }
}

/// Calcula nutação em longitude (`nutl`) e em obliquidade (`nuto`) em radianos para a data `jd` (`nutate.c`).
pub fn nutlo(jd: f64) -> (f64, f64) {
    let t = (jd - J2000) / 36525.0;
    let t2 = t * t;
    let t10 = t / 10.0;

    let om = (mod3600(-6962890.539 * t + 450160.280) + (0.008 * t + 7.455) * t2) * STR;
    let ms = (mod3600(129596581.224 * t + 1287099.804) - (0.012 * t + 0.577) * t2) * STR;
    let mm = (mod3600(1717915922.633 * t + 485866.733) + (0.064 * t + 31.310) * t2) * STR;
    let ff = (mod3600(1739527263.137 * t + 335778.877) + (0.011 * t - 13.257) * t2) * STR;
    let dd = (mod3600(1602961601.328 * t + 1072261.307) + (0.019 * t - 6.891) * t2) * STR;

    let mut ss = [[0.0f64; 8]; 5];
    let mut cc = [[0.0f64; 8]; 5];
    fill_sscc(&mut ss, &mut cc, 0, mm, 3);
    fill_sscc(&mut ss, &mut cc, 1, ms, 2);
    fill_sscc(&mut ss, &mut cc, 2, ff, 4);
    fill_sscc(&mut ss, &mut cc, 3, dd, 4);
    fill_sscc(&mut ss, &mut cc, 4, om, 2);

    let mut c_sum = 0.0;
    let mut d_sum = 0.0;

    for i in 0..105 {
        let base = i * 9;
        let mut k1 = false;
        let mut cv = 0.0;
        let mut sv = 0.0;

        for m in 0..5 {
            let j = NT[base + m];
            if j != 0 {
                let k = j.unsigned_abs() as usize;
                let mut su = ss[m][k - 1];
                if j < 0 {
                    su = -su;
                }
                let cu = cc[m][k - 1];
                if !k1 {
                    sv = su;
                    cv = cu;
                    k1 = true;
                } else {
                    let sw = su * cv + cu * sv;
                    cv = cu * cv - su * sv;
                    sv = sw;
                }
            }
        }

        let mut f = NT[base + 5] as f64;
        let k_lon = NT[base + 6];
        if k_lon != 0 {
            f += t10 * (k_lon as f64);
        }

        let mut g = NT[base + 7] as f64;
        let k_obl = NT[base + 8];
        if k_obl != 0 {
            g += t10 * (k_obl as f64);
        }

        c_sum += f * sv;
        d_sum += g * cv;
    }

    c_sum += (-1742.0 * t10 - 171996.0) * ss[4][0];
    d_sum += (89.0 * t10 + 92025.0) * cc[4][0];

    let nutl = 0.0001 * STR * c_sum;
    let nuto = 0.0001 * STR * d_sum;
    (nutl, nuto)
}

/// Aplica rotação de nutação ao vetor equatorial retangular `p` na data `jd` (`nutate.c`).
pub fn nutate(jd: f64, p: &mut [f64; 3]) {
    let (nutl, nuto) = nutlo(jd);
    let (eps, coseps, sineps) = epsiln(jd);

    let f = eps + nuto;
    let ce = f.cos();
    let se = f.sin();
    let sino = nuto.sin();
    let cl = nutl.cos();
    let sl = nutl.sin();

    let p0 = cl * p[0] - sl * coseps * p[1] - sl * sineps * p[2];
    let p1 = sl * ce * p[0]
        + (cl * coseps * ce + sineps * se) * p[1]
        - (sino + (1.0 - cl) * sineps * ce) * p[2];
    let p2 = sl * se * p[0]
        + (sino + (cl - 1.0) * se * coseps) * p[1]
        + (cl * sineps * se + coseps * ce) * p[2];

    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
}

// Coeficientes de precessão de Williams (1994) introduzidos no DE403 (`precess.c` linhas 85-86 do aa-56):
const PA_COF: [f64; 10] = [
    -8.66e-10, -4.759e-8, 2.424e-7, 1.3095e-5, 1.7451e-4, -1.8055e-3, -0.235316, 0.076, 110.5414,
    50287.91959,
];

const NODE_COF: [f64; 11] = [
    6.6402e-16,
    -2.69151e-15,
    -1.547021e-12,
    7.521313e-12,
    1.9e-10,
    -3.54e-9,
    -1.8103e-7,
    1.26e-7,
    7.436169e-5,
    -0.04207794833,
    3.052115282424,
];

const INCL_COF: [f64; 11] = [
    1.2147e-16,
    7.3759e-17,
    -8.26287e-14,
    2.503410e-13,
    2.4650839e-11,
    -5.4000441e-11,
    1.32115526e-9,
    -6.012e-7,
    -1.62442e-5,
    0.00227850649,
    0.0,
];

/// Precessão do equinócio e eclíptica (`precess.c` - método Williams DE403 do `aa-56`).
/// `direction = -1`: de J2000 para `jd`; `direction = 1`: de `jd` para J2000.
pub fn precess(r: &mut [f64; 3], jd: f64, direction: i32) {
    if (jd - J2000).abs() < 1e-12 {
        return;
    }

    let (_, coseps_init, sineps_init) = if direction == 1 {
        epsiln(jd)
    } else {
        epsiln(J2000)
    };

    let t = (jd - J2000) / 365250.0; // milhares de anos Julianos

    let mut pa = PA_COF[0];
    for &c in &PA_COF[1..] {
        pa = pa * t + c;
    }
    pa *= STR * t;

    let mut w = NODE_COF[0];
    for &c in &NODE_COF[1..] {
        w = w * t + c;
    }

    let mut th = INCL_COF[0];
    for &c in &INCL_COF[1..] {
        th = th * t + c;
    }

    let mut x = [0.0; 3];
    x[0] = r[0];
    let z_tmp = coseps_init * r[1] + sineps_init * r[2];
    x[2] = -sineps_init * r[1] + coseps_init * r[2];
    x[1] = z_tmp;

    let z1 = if direction == 1 { w + pa } else { w };
    let b1 = z1.cos();
    let a1 = z1.sin();
    let z_tmp = b1 * x[0] + a1 * x[1];
    x[1] = -a1 * x[0] + b1 * x[1];
    x[0] = z_tmp;

    let th_rot = if direction == 1 { -th } else { th };
    let b2 = th_rot.cos();
    let a2 = th_rot.sin();
    let z_tmp = b2 * x[1] + a2 * x[2];
    x[2] = -a2 * x[1] + b2 * x[2];
    x[1] = z_tmp;

    let z3 = if direction == 1 { -w } else { -w - pa };
    let b3 = z3.cos();
    let a3 = z3.sin();
    let z_tmp = b3 * x[0] + a3 * x[1];
    x[1] = -a3 * x[0] + b3 * x[1];
    x[0] = z_tmp;

    let (_, coseps_fin, sineps_fin) = if direction == 1 {
        epsiln(J2000)
    } else {
        epsiln(jd)
    };
    let z_tmp = coseps_fin * x[1] - sineps_fin * x[2];
    x[2] = sineps_fin * x[1] + coseps_fin * x[2];
    x[1] = z_tmp;

    *r = x;
}

/// Converte vetor retangular equatorial `pp` em coordenadas polares eclípticas `(lon, lat, r)` (`lonlat.c`).
/// Se `ofdate` for `true` (`1` em C), precessiona `pp` de J2000 para `jd`; caso contrário, assume que `pp` já está em `jd` e usa `epsiln(jd)`.
pub fn lonlat(pp: &[f64; 3], jd: f64, ofdate: bool) -> (f64, f64, f64) {
    let mut s = *pp;
    let r = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
    if ofdate {
        precess(&mut s, jd, -1);
    }
    let (_, coseps, sineps) = epsiln(jd);

    let yy = s[1];
    let zz = s[2];
    let x = s[0];
    let y = coseps * yy + sineps * zz;
    let z = -sineps * yy + coseps * zz;

    let lon = zatan2(x, y);
    let lat = if r > 0.0 {
        (z / r).clamp(-1.0, 1.0).asin()
    } else {
        0.0
    };
    (lon, lat, r)
}

/// Tempo Sideral Aparente Local em segundos (`sidrlt.c` do `aa-56` - fórmula Williams 1994 / DE403).
pub fn sidrlt(jd_ut: f64, tdt: f64, tlong_deg: f64) -> f64 {
    let mut jd0 = jd_ut.floor();
    let mut secs = jd_ut - jd0;
    if secs < 0.5 {
        jd0 -= 0.5;
        secs += 0.5;
    } else {
        jd0 += 0.5;
        secs -= 0.5;
    }
    secs *= 86400.0;

    let t0 = (jd0 - J2000) / 36525.0;
    let (nutl, _) = nutlo(tdt);
    let (_, coseps, _) = epsiln(tdt);
    let eqeq = 240.0 * RTD * nutl * coseps;

    // Corrections to Williams (1994) introduced in DE403 (sidrlt.c lines 78-81):
    let gmst0 = (((-2.0e-6 * t0 - 3.0e-7) * t0 + 9.27701e-2) * t0 + 8640184.7942063) * t0
        + 24110.54841;
    let msday = (((-(4.0 * 2.0e-6) * t0 - (3.0 * 3.0e-7)) * t0 + (2.0 * 9.27701e-2)) * t0
        + 8640184.7942063)
        / (86400.0 * 36525.0)
        + 1.0;

    let gmst = gmst0 + msday * secs + eqeq + 240.0 * tlong_deg;
    gmst - 86400.0 * (gmst / 86400.0).floor()
}

/// Refração atmosférica em graus (`refrac.c`).
pub fn refrac(alt_deg: f64, atpress: f64, attemp: f64) -> f64 {
    if alt_deg < -2.0 || alt_deg >= 90.0 {
        return 0.0;
    }
    if alt_deg > 15.0 {
        return 0.00452 * atpress / ((273.0 + attemp) * (DTR * alt_deg).tan());
    }
    let mut y = alt_deg;
    let mut d = 0.0;
    let p = (atpress - 80.0) / 930.0;
    let q = 4.8e-3 * (attemp - 10.0);
    let mut y0 = y;
    let mut d0 = d;

    for _ in 0..4 {
        let mut n = y + (7.31 / (y + 4.4));
        n = 1.0 / (DTR * n).tan();
        d = n * p / (60.0 + q * (n + 39.0));
        n = y - y0;
        y0 = d - d0 - n;
        let next_y = if n != 0.0 && y0 != 0.0 {
            y - n * (alt_deg + d - y) / y0
        } else {
            alt_deg + d
        };
        y0 = y;
        d0 = d;
        y = next_y;
    }
    d
}

/// Redução de latitude geodésica `glat_deg` e altura `height_m` para latitude geocêntrica `tlat_deg` e raio terrestre `trho` (`kfiles.c` e `moonrise.c`).
pub fn geodetic_to_geocentric(glat_deg: f64, height_m: f64) -> (f64, f64) {
    let u = glat_deg * DTR;
    let co = u.cos();
    let si = u.sin();
    let fl = (1.0 - 1.0 / FLAT_INV).powi(2);
    let si2 = si * si;
    let u_inv = 1.0 / (co * co + fl * si2).sqrt();
    let a = AEARTH_M * u_inv + height_m;
    let b = AEARTH_M * fl * u_inv + height_m;
    let trho_m = (a * a * co * co + b * b * si2).sqrt();
    let mut tlat = RTD * (a * co / trho_m).clamp(-1.0, 1.0).acos();
    if glat_deg < 0.0 {
        tlat = -tlat;
    }
    let trho = trho_m / AEARTH_M;
    (tlat, trho)
}

/// Calcula Altitude e Azimute Topocêntricos (`altaz.c`, `diurab.c`, `diurpx.c`).
pub fn topocentric_altaz(
    ra_rad: f64,
    dec_rad: f64,
    dist_au: f64,
    jd_ut: f64,
    tdt: f64,
    obs_lon_deg: f64,
    obs_lat_deg: f64,
) -> (f64, f64, f64) {
    let last_sec = sidrlt(jd_ut, tdt, obs_lon_deg);
    let last = last_sec * DTR / 240.0;

    let (tlat_deg, trho) = geodetic_to_geocentric(obs_lat_deg, 0.0);
    let tlat_rad = tlat_deg * DTR;

    let mut ra = ra_rad;
    let mut dec = dec_rad;

    // Diurnal aberration (`diurab.c`)
    let lha = last - ra;
    let cosdec = dec.cos();
    let sindec = dec.sin();
    let coslat_gc = tlat_rad.cos();
    let sinlat_gc = tlat_rad.sin();
    if cosdec.abs() > 1e-12 {
        ra += 1.5472e-6 * trho * coslat_gc * lha.cos() / cosdec;
    }
    dec += 1.5472e-6 * trho * coslat_gc * lha.sin() * sindec;

    // Diurnal parallax (`diurpx.c`)
    if dist_au < 1758.8 && dist_au > 0.0 {
        let disfac = AU_KM / (0.001 * AEARTH_M);
        let d = dist_au * disfac;
        let x = d * dec.cos() * ra.cos() - trho * coslat_gc * last.cos();
        let y = d * dec.cos() * ra.sin() - trho * coslat_gc * last.sin();
        let z = d * dec.sin() - trho * sinlat_gc;
        let top_d = (x * x + y * y + z * z).sqrt();
        ra = zatan2(x, y);
        dec = (z / top_d).clamp(-1.0, 1.0).asin();
    }

    let lha = last - ra;
    let cosdec = dec.cos();
    let sindec = dec.sin();
    let coslha = lha.cos();
    let sinlha = lha.sin();
    let glat_rad = DTR * obs_lat_deg;
    let coslat = glat_rad.cos();
    let sinlat = glat_rad.sin();

    let n = -cosdec * sinlha;
    let d_az = sindec * coslat - cosdec * coslha * sinlat;
    let az_deg = RTD * zatan2(d_az, n);
    let alt_true_deg = RTD * (sindec * sinlat + cosdec * coslha * coslat).clamp(-1.0, 1.0).asin();
    let refr = refrac(alt_true_deg, 1013.0, 20.0);
    let alt_app_deg = alt_true_deg + refr;

    (az_deg, alt_true_deg, alt_app_deg)
}

