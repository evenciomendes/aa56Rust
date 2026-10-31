//! CLI entry point (`aa56-cli`) for `aa56-rust`.
//! Supports JSON output for `--mode ephemeris | stars | stars1950 | messier | orbits | conjunct | moonrise`.

use serde::Serialize;
use std::env;

use aa56_rust::conjunct::compute_conjunct_report;
use aa56_rust::kepler::{kepler, De404Body, OrbitElement};
use aa56_rust::math_time::{caltoj, deltat};
use aa56_rust::moonrise::compute_moonrise_table;
use aa56_rust::reductions::{compute_ephemeris_report, EphemerisInput};
use aa56_rust::stars_catalog::{
    parse_star_catalog, reduce_orbit_catalog, reduce_star, ReducedOrbitEntry, ReducedStarEntry,
};
use aa56_rust::{CATALOG_MESSIER, CATALOG_ORBIT, CATALOG_STAR, CATALOG_STAR_1950};

#[derive(Serialize)]
struct StarsCatalogResponse {
    catalog: String,
    count: usize,
    jd_ut: f64,
    jd_tdt: f64,
    delta_t_sec: f64,
    items: Vec<ReducedStarEntry>,
}

#[derive(Serialize)]
struct OrbitsCatalogResponse {
    catalog: String,
    count: usize,
    jd_ut: f64,
    jd_tdt: f64,
    delta_t_sec: f64,
    items: Vec<ReducedOrbitEntry>,
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let mut mode = "ephemeris".to_string();
    let mut year: i32 = 2025;
    let mut month: i32 = 1;
    let mut day: f64 = 1.0;
    let mut hour: f64 = 12.0;
    let mut minute: f64 = 0.0;
    let mut second: f64 = 0.0;
    let mut dt_override: Option<f64> = None;
    let mut lon: f64 = -71.13;
    let mut lat: f64 = 42.38;
    let mut height: f64 = 24.0;
    let mut temp: f64 = 12.0;
    let mut press: f64 = 1010.0;
    let mut days_span: usize = 14;
    let mut cycles_span: usize = 3;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--mode" if i + 1 < args.len() => {
                mode = args[i + 1].clone();
                i += 2;
            }
            "--year" if i + 1 < args.len() => {
                year = args[i + 1].parse().unwrap_or(2025);
                i += 2;
            }
            "--month" if i + 1 < args.len() => {
                month = args[i + 1].parse().unwrap_or(1);
                i += 2;
            }
            "--day" if i + 1 < args.len() => {
                day = args[i + 1].parse().unwrap_or(1.0);
                i += 2;
            }
            "--hour" if i + 1 < args.len() => {
                hour = args[i + 1].parse().unwrap_or(0.0);
                i += 2;
            }
            "--minute" if i + 1 < args.len() => {
                minute = args[i + 1].parse().unwrap_or(0.0);
                i += 2;
            }
            "--second" if i + 1 < args.len() => {
                second = args[i + 1].parse().unwrap_or(0.0);
                i += 2;
            }
            "--dt" if i + 1 < args.len() => {
                if let Ok(v) = args[i + 1].parse::<f64>() {
                    dt_override = Some(v);
                }
                i += 2;
            }
            "--lon" if i + 1 < args.len() => {
                lon = args[i + 1].parse().unwrap_or(-71.13);
                i += 2;
            }
            "--lat" if i + 1 < args.len() => {
                lat = args[i + 1].parse().unwrap_or(42.38);
                i += 2;
            }
            "--height" if i + 1 < args.len() => {
                height = args[i + 1].parse().unwrap_or(24.0);
                i += 2;
            }
            "--temp" if i + 1 < args.len() => {
                temp = args[i + 1].parse().unwrap_or(12.0);
                i += 2;
            }
            "--press" if i + 1 < args.len() => {
                press = args[i + 1].parse().unwrap_or(1010.0);
                i += 2;
            }
            "--days" if i + 1 < args.len() => {
                days_span = args[i + 1].parse().unwrap_or(14);
                i += 2;
            }
            "--cycles" if i + 1 < args.len() => {
                cycles_span = args[i + 1].parse().unwrap_or(3);
                i += 2;
            }
            _ => {
                i += 1;
            }
        }
    }

    let day_frac = day + (hour + minute / 60.0 + second / 3600.0) / 24.0;
    let jd_ut = caltoj(year, month, day_frac);
    let y_dec = (year as f64) + ((month as f64) - 1.0) / 12.0 + day_frac / 365.25;
    let delta_t_sec = dt_override.unwrap_or_else(|| deltat(y_dec));
    let jd_tdt = jd_ut + delta_t_sec / 86400.0;

    match mode.as_str() {
        "ephemeris" => {
            let input = EphemerisInput {
                year,
                month,
                day,
                hour,
                minute,
                second,
                delta_t_override: dt_override,
                obs_lon_deg: lon,
                obs_lat_deg: lat,
                obs_height_m: height,
                temp_c: temp,
                pressure_mb: press,
            };
            let report = compute_ephemeris_report(&input);
            println!("{}", serde_json::to_string_pretty(&report).unwrap());
        }
        "stars" | "stars1950" | "messier" => {
            let (cat_name, raw_cat) = match mode.as_str() {
                "stars1950" => ("star1950.cat (FK4 B1950 -> FK5 J2000)", CATALOG_STAR_1950),
                "messier" => ("messier.cat (Messier Deep-Sky Objects)", CATALOG_MESSIER),
                _ => ("star.cat (aa-56 Navigation & Reference Stars)", CATALOG_STAR),
            };
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
            let (earth_j2000, _) = kepler(jd_tdt, &mut earth_el);
            let recs = parse_star_catalog(raw_cat);
            let items: Vec<ReducedStarEntry> = recs
                .iter()
                .enumerate()
                .map(|(idx, r)| {
                    reduce_star(
                        idx + 1,
                        r,
                        &earth_j2000,
                        jd_tdt,
                        jd_ut,
                        lon,
                        lat,
                        height,
                        temp,
                        press,
                    )
                })
                .collect();
            let resp = StarsCatalogResponse {
                catalog: cat_name.to_string(),
                count: items.len(),
                jd_ut,
                jd_tdt,
                delta_t_sec,
                items,
            };
            println!("{}", serde_json::to_string_pretty(&resp).unwrap());
        }
        "orbits" => {
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
            let (earth_j2000, _) = kepler(jd_tdt, &mut earth_el);
            let items = reduce_orbit_catalog(
                CATALOG_ORBIT,
                &earth_j2000,
                jd_tdt,
                jd_ut,
                lon,
                lat,
                height,
                temp,
                press,
            );
            let resp = OrbitsCatalogResponse {
                catalog: "orbit.cat (aa-56 Comets & Minor Planets)".to_string(),
                count: items.len(),
                jd_ut,
                jd_tdt,
                delta_t_sec,
                items,
            };
            println!("{}", serde_json::to_string_pretty(&resp).unwrap());
        }
        "conjunct" => {
            let rep = compute_conjunct_report(jd_tdt, cycles_span);
            println!("{}", serde_json::to_string_pretty(&rep).unwrap());
        }
        "moonrise" => {
            let rep = compute_moonrise_table(jd_ut, days_span, lon, lat, height, temp, press);
            println!("{}", serde_json::to_string_pretty(&rep).unwrap());
        }
        _ => {
            eprintln!("Unknown mode: {}", mode);
            std::process::exit(1);
        }
    }
}

