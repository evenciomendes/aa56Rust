//! Constantes fundamentais e utilitários matemáticos portados 100% de `aa.c`, `kfiles.c`, `dms.c`, `zatan2.c` e `deltat.c` do `aa-56`.

use std::f64::consts::PI;

pub const DTR: f64 = 1.7453292519943295769e-2; // graus para radianos (aa.c)
pub const RTD: f64 = 5.7295779513082320877e1;  // radianos para graus (aa.c)
pub const RTS: f64 = 2.0626480624709635516e5;  // segundos de arco por radiano (aa.c)
pub const STR: f64 = 4.8481368110953599359e-6; // radianos por segundo de arco (aa.c)
pub const J2000: f64 = 2451545.0;              // 2000 January 1.5 (aa.c)
pub const B1950: f64 = 2433282.423;            // 1950 January 0.923 Besselian epoch (aa.c)
pub const J1900: f64 = 2415020.0;              // 1900 January 0, 12h UT (aa.c)

// Constantes DE403/DE404 definidas em kfiles.c:
pub const AEARTH_M: f64 = 6378137.0;           // Radius of the earth, in meters (kfiles.c)
pub const AU_KM: f64 = 1.49597870691e8;        // Astronomical unit, in kilometers (kfiles.c)
pub const EMRAT: f64 = 81.300585;              // Earth/Moon mass ratio (kfiles.c)
pub const CLIGHT_KMS: f64 = 2.99792458e5;      // Speed of light, km/sec (kfiles.c)
pub const CLIGHT_AUD: f64 = 86400.0 * CLIGHT_KMS / AU_KM; // ~173.1446326846693 AU/day (kfiles.c)
pub const REARTH_AU: f64 = 0.001 * AEARTH_M / AU_KM;      // ~4.263523250891372e-5 AU (kfiles.c)
pub const FLAT_INV: f64 = 298.257222;          // Flattening factor (kfiles.c)

/// Equivalente exato de `zatan2(x, y)` de Stephen Moshier (`zatan2.c`).
/// Atenção: No aa-56, `zatan2(x, y)` computa `atan(y / x)` no intervalo `[0, 2*PI)`.
pub fn zatan2(x: f64, y: f64) -> f64 {
    let mut code: u8 = 0;
    if x < 0.0 {
        code = 2;
    }
    if y < 0.0 {
        code |= 1;
    }
    if x == 0.0 {
        if (code & 1) != 0 {
            return 1.5 * PI;
        }
        if y == 0.0 {
            return 0.0;
        }
        return 0.5 * PI;
    }
    if y == 0.0 {
        if (code & 2) != 0 {
            return PI;
        }
        return 0.0;
    }
    let w = match code {
        0 => 0.0,
        1 => 2.0 * PI,
        2 | 3 => PI,
        _ => 0.0,
    };
    w + (y / x).atan()
}

pub fn mod360(x: f64) -> f64 {
    let k = (x / 360.0) as i64;
    let mut y = x - (k as f64) * 360.0;
    while y < 0.0 {
        y += 360.0;
    }
    while y > 360.0 {
        y -= 360.0;
    }
    y
}

pub fn modtp(x: f64) -> f64 {
    let tpi = 2.0 * PI;
    let mut y = (x / tpi).floor();
    y = x - y * tpi;
    while y < 0.0 {
        y += tpi;
    }
    while y >= tpi {
        y -= tpi;
    }
    y
}

pub fn mod3600(x: f64) -> f64 {
    x - 1296000.0 * (x / 1296000.0).floor()
}

pub fn mods3600(x: f64) -> f64 {
    x - 1.296e6 * (x / 1.296e6).floor()
}

/// Converte mudança em coordenadas retangulares para mudança em RA e Dec (`deltap` em `dms.c`).
pub fn deltap(p0: &[f64; 3], p1: &[f64; 3]) -> (f64, f64) {
    let mut dp = [0.0; 3];
    let mut p_sq = 0.0;
    let mut q_sq = 0.0;
    let mut z = 0.0;
    for i in 0..3 {
        let x = p0[i];
        let y = p1[i];
        p_sq += x * x;
        q_sq += y * y;
        let diff = y - x;
        dp[i] = diff;
        z += diff * diff;
    }
    let a = p_sq.sqrt();
    let b = q_sq.sqrt();
    if a < 1.0e-7 || b < 1.0e-7 || (z / (p_sq + q_sq)) > 5.0e-7 {
        let p_ra = zatan2(p0[0], p0[1]);
        let mut q_ra = zatan2(p1[0], p1[1]) - p_ra;
        while q_ra < -PI {
            q_ra += 2.0 * PI;
        }
        while q_ra > PI {
            q_ra -= 2.0 * PI;
        }
        let p_dec = (p0[2] / a).clamp(-1.0, 1.0).asin();
        let q_dec = (p1[2] / b).clamp(-1.0, 1.0).asin();
        return (q_ra, q_dec - p_dec);
    }
    let x = p0[0];
    let y = p0[1];
    let dr = if x == 0.0 {
        1.0e38
    } else {
        let q = y / x;
        (dp[1] - dp[0] * y / x) / (x * (1.0 + q * q))
    };
    let x_dec = p0[2] / a;
    let p_cos = (1.0 - x_dec * x_dec).max(1e-18).sqrt();
    let dd = (p1[2] / b - x_dec) / p_cos;
    (dr, dd)
}

/// Converte data calendário Gregoriano/Juliano para Dia Juliano (`caltoj` em `dms.c`).
pub fn caltoj(year: impl Into<i64>, month: i32, day: f64) -> f64 {
    let year = year.into();
    let mut y = year + 4800;
    if year < 0 {
        y += 1;
    }
    let mut m = month as i64;
    if m <= 2 {
        m += 12;
        y -= 1;
    }
    let e = (306 * (m + 1)) / 10;
    let a = y / 100;

    let is_julian = if year < 1582 {
        true
    } else if year == 1582 {
        if month < 10 {
            true
        } else if month > 10 {
            false
        } else {
            day < 15.0
        }
    } else {
        false
    };

    let b = if is_julian {
        -38
    } else {
        (a / 4) - a
    };

    let c = (36525 * y) / 100;
    (b + c + e) as f64 + day - 32167.5
}

/// Converte Dia Juliano para data calendário Gregoriano/Juliano (`jtocal` em `dms.c`).
/// Retorna `(year, month, day_int, hour, minute, second, formatted_string)`.
pub fn jtocal_parts(j: f64) -> (i64, i32, i32, i32, i32, f64, String) {
    let bc = j < 1721423.5;
    let jd = (j + 0.5).floor() as i64;
    let a = (100 * jd + 3204500) / 3652425;
    let mut c = jd + 1486;
    if (jd as f64) >= 2299160.5 {
        c += a - a / 4;
    } else {
        c += 38;
    }
    let d = (100 * c - 12210) / 36525;
    let x = (36525 * d) / 100;
    let y = ((c - x) * 100) / 3061;
    let day_calc = (c - x - ((306 * y) / 10)) as i32;
    let mut month = (y - 1) as i32;
    if y > 13 {
        month -= 12;
    }
    let mut year = d - 4715;
    if month > 2 {
        year -= 1;
    }
    let dd = (day_calc as f64) + j - (jd as f64) + 0.5;
    let day_int = dd.floor() as i32;
    let frac = (dd - (day_int as f64)).clamp(0.0, 0.999999999);
    let total_sec = frac * 86400.0;
    let h = (total_sec / 3600.0).floor() as i32;
    let rem = total_sec - (h as f64) * 3600.0;
    let m = (rem / 60.0).floor() as i32;
    let s = rem - (m as f64) * 60.0;

    let disp_year = if bc { -year + 1 } else { year };
    let str_date = if bc {
        format!("{:04} B.C.-{:02}-{:02} {:02}:{:02}:{:04.1}", disp_year, month, day_int, h, m, s)
    } else {
        format!("{:04}-{:02}-{:02} {:02}:{:02}:{:04.1}", disp_year, month, day_int, h, m, s)
    };
    (year, month, day_int, h, m, s, str_date)
}

/// Formata graus em `GG° MM' SS.ss″` (`dms` em `dms.c`)
pub fn format_dms(deg: f64) -> String {
    let sign = if deg < 0.0 { "-" } else { "" };
    let mut s = deg.abs();
    let d = s.floor() as i32;
    s = (s - d as f64) * 60.0;
    let m = s.floor() as i32;
    s = (s - m as f64) * 60.0;
    format!("{}{:02}° {:02}' {:05.2}″", sign, d, m, s)
}

/// Formata radianos em horas, minutos e segundos `HHh MMm SS.sss` (`hms` em `dms.c`)
pub fn format_hms_rad(x: f64) -> String {
    let rtoh = 12.0 / PI;
    let mut s = x * rtoh;
    while s < 0.0 {
        s += 24.0;
    }
    while s >= 24.0 {
        s -= 24.0;
    }
    let mut h = s.floor() as i32;
    s = (s - h as f64) * 60.0;
    let mut m = s.floor() as i32;
    s = (s - m as f64) * 60.0;
    let mut sfrac = (1000.0 * s + 0.5) as i64;
    if sfrac >= 60000 {
        sfrac -= 60000;
        m += 1;
        if m >= 60 {
            m -= 60;
            h = (h + 1) % 24;
        }
    }
    let sint = sfrac / 1000;
    let rem = sfrac - sint * 1000;
    format!("{:02}h {:02}m {:02}.{:03}s", h, m, sint, rem)
}

/// Tabela Morrison & Stephenson (2004) de -1000 a 1700 em passos de 100 anos (`deltat.c` do `aa-56`)
const MS_TABLE: [i16; 28] = [
    25428, 23700, 22000, 21000, 19040, 17190, 15530, 14080, 12790, 11640,
    10580, 9600, 8640, 7680, 6700, 5710, 4740, 3810, 2960, 2200,
    1570, 1090, 740, 490, 320, 200, 120, 9,
];

/// Tabela de DeltaT de 1620.0 a 2017.0 (398 entradas) em centésimos de segundo (`deltat.c` do `aa-56`)
const DT_TABLE: [i16; 398] = [
    12400, 11900, 11500, 11000, 10600, 10200, 9800, 9500, 9100, 8800,
    8500, 8200, 7900, 7700, 7400, 7200, 7000, 6700, 6500, 6300,
    6200, 6000, 5800, 5700, 5500, 5400, 5300, 5100, 5000, 4900,
    4800, 4700, 4600, 4500, 4400, 4300, 4200, 4100, 4000, 3800,
    3700, 3600, 3500, 3400, 3300, 3200, 3100, 3000, 2800, 2700,
    2600, 2500, 2400, 2300, 2200, 2100, 2000, 1900, 1800, 1700,
    1600, 1500, 1400, 1400, 1300, 1200, 1200, 1100, 1100, 1000,
    1000, 1000, 900, 900, 900, 900, 900, 900, 900, 900,
    900, 900, 900, 900, 900, 900, 900, 900, 1000, 1000,
    1000, 1000, 1000, 1000, 1000, 1000, 1000, 1100, 1100, 1100,
    1100, 1100, 1100, 1100, 1100, 1100, 1100, 1100, 1100, 1100,
    1100, 1100, 1100, 1100, 1200, 1200, 1200, 1200, 1200, 1200,
    1200, 1200, 1200, 1200, 1300, 1300, 1300, 1300, 1300, 1300,
    1300, 1400, 1400, 1400, 1400, 1400, 1400, 1400, 1500, 1500,
    1500, 1500, 1500, 1500, 1500, 1600, 1600, 1600, 1600, 1600,
    1600, 1600, 1600, 1600, 1600, 1700, 1700, 1700, 1700, 1700,
    1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
    1700, 1700, 1600, 1600, 1600, 1600, 1500, 1500, 1400, 1400,
    1370, 1340, 1310, 1290, 1270, 1260, 1250, 1250, 1250, 1250,
    1250, 1250, 1250, 1250, 1250, 1250, 1250, 1240, 1230, 1220,
    1200, 1170, 1140, 1110, 1060, 1020, 960, 910, 860, 800,
    750, 700, 660, 630, 600, 580, 570, 560, 560, 560,
    570, 580, 590, 610, 620, 630, 650, 660, 680, 690,
    710, 720, 730, 740, 750, 760, 770, 770, 780, 780,
    788, 782, 754, 697, 640, 602, 541, 410, 292, 182,
    161, 10, -102, -128, -269, -324, -364, -454, -471, -511,
    -540, -542, -520, -546, -546, -579, -563, -564, -580, -566,
    -587, -601, -619, -664, -644, -647, -609, -576, -466, -374,
    -272, -154, -2, 124, 264, 386, 537, 614, 775, 913,
    1046, 1153, 1336, 1465, 1601, 1720, 1824, 1906, 2025, 2095,
    2116, 2225, 2241, 2303, 2349, 2362, 2386, 2449, 2434, 2408,
    2402, 2400, 2387, 2395, 2386, 2393, 2373, 2392, 2396, 2402,
    2433, 2483, 2530, 2570, 2624, 2677, 2728, 2778, 2825, 2871,
    2915, 2957, 2997, 3036, 3072, 3107, 3135, 3168, 3218, 3268,
    3315, 3359, 3400, 3447, 3503, 3573, 3654, 3743, 3829, 3920,
    4018, 4117, 4223, 4337, 4449, 4548, 4646, 4752, 4853, 4959,
    5054, 5138, 5217, 5296, 5379, 5434, 5487, 5532, 5582, 5630,
    5686, 5757, 5831, 5912, 5998, 6078, 6163, 6230, 6297, 6347,
    6383, 6409, 6430, 6447, 6457, 6469, 6485, 6515, 6546, 6578,
    6607, 6632, 6660, 6691, 6728, 6764, 6810, 6859,
];

/// Calcula `DeltaT = TDT - UT` em segundos para um ano Juliano `y` (`deltat.c` do `aa-56`, `TABEND = 2017`).
pub fn deltat(y: f64) -> f64 {
    const TABSTART: i32 = 1620;
    const TABEND: f64 = 2017.0;
    const TABSIZ: usize = 398;

    if y > TABEND {
        if y > TABEND + 100.0 {
            let b = 0.01 * (y - 1820.0);
            return 32.0 * b * b - 20.0;
        } else {
            let a = 0.01 * (DT_TABLE[TABSIZ - 1] as f64);
            let b = 0.001 * ((DT_TABLE[TABSIZ - 1] - DT_TABLE[TABSIZ - 11]) as f64);
            let b_lt = 0.01 * (TABEND + 100.0 - 1820.0);
            let m0 = 32.0 * b_lt * b_lt - 20.0;
            let m1 = 0.64 * b_lt;
            let d = 2.0e-6 * (50.0 * (m1 + b) - m0 + a);
            let c = 1.0e-4 * (m0 - a - 100.0 * b - 1.0e6 * d);
            let p = y - TABEND;
            return a + p * (b + p * (c + p * d));
        }
    }

    if y < 1700.0 {
        if y <= -1000.0 {
            let b = 0.01 * (y - 1820.0);
            return 32.0 * b * b - 20.0;
        } else {
            let iy = ((y as i32) + 1000) / 100;
            let idx = iy.clamp(0, 26) as usize;
            let b = -1000.0 + 100.0 * (idx as f64);
            let p = MS_TABLE[idx] as f64;
            return p + 0.01 * (y - b) * ((MS_TABLE[idx + 1] as f64) - p);
        }
    }

    let p_floor = y.floor();
    let iy = (p_floor as i32) - TABSTART;
    if iy < 0 {
        return 0.01 * (DT_TABLE[0] as f64);
    }
    let iy_u = iy as usize;
    let mut ans = DT_TABLE[iy_u] as f64;
    let k = iy_u + 1;
    if k >= TABSIZ {
        return 0.01 * ans;
    }

    let p = y - p_floor;
    ans += p * ((DT_TABLE[k] - DT_TABLE[iy_u]) as f64);
    if iy - 1 < 0 || iy_u + 2 >= TABSIZ {
        return 0.01 * ans;
    }

    let mut d = [0i32; 6];
    let mut kd = iy - 2;
    for i in 0..5 {
        if kd < 0 || (kd + 1) >= (TABSIZ as i32) {
            d[i] = 0;
        } else {
            d[i] = (DT_TABLE[(kd + 1) as usize] - DT_TABLE[kd as usize]) as i32;
        }
        kd += 1;
    }
    for i in 0..4 {
        d[i] = d[i + 1] - d[i];
    }
    let mut b = 0.25 * p * (p - 1.0);
    ans += b * ((d[1] + d[2]) as f64);
    if iy_u + 2 >= TABSIZ {
        return 0.01 * ans;
    }
    for i in 0..3 {
        d[i] = d[i + 1] - d[i];
    }
    b = 2.0 * b / 3.0;
    ans += (p - 0.5) * b * (d[1] as f64);
    if iy - 2 < 0 || iy_u + 3 > TABSIZ {
        return 0.01 * ans;
    }
    for i in 0..2 {
        d[i] = d[i + 1] - d[i];
    }
    b = 0.125 * b * (p + 1.0) * (p - 2.0);
    ans += b * ((d[0] + d[1]) as f64);

    0.01 * ans
}
