//! Daily / Multi-latitude Moonrise, Transit, and Moonset Tables (`moonrise.c`).
//! Ports Stephen L. Moshier's `moonrise.c` to generate daily lunar rise, culmination, and set tables.

use serde::Serialize;

use crate::domoon::domoon;
use crate::math_time::{deltat, jtocal_parts};

#[derive(Debug, Clone, Serialize)]
pub struct MoonriseDayEntry {
    pub day_index: usize,
    pub date_str: String,
    pub jd_ut_noon: f64,
    pub rise_ut: Option<String>,
    pub transit_ut: Option<String>,
    pub set_ut: Option<String>,
    pub rise_jd: Option<f64>,
    pub transit_jd: Option<f64>,
    pub set_jd: Option<f64>,
    pub illuminated_pct: f64,
    pub phase_description: String,
    pub app_ra_hms: String,
    pub app_dec_dms: String,
    pub horiz_parallax_dms: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct MoonriseReport {
    pub start_jd_ut: f64,
    pub num_days: usize,
    pub obs_lon_deg: f64,
    pub obs_lat_deg: f64,
    pub obs_height_m: f64,
    pub entries: Vec<MoonriseDayEntry>,
}

/// Compute a multi-day Moonrise / Transit / Moonset table (`moonrise.c`).
pub fn compute_moonrise_table(
    start_jd_ut: f64,
    num_days: usize,
    obs_lon_deg: f64,
    obs_lat_deg: f64,
    obs_height_m: f64,
    temp_c: f64,
    pressure_mb: f64,
) -> MoonriseReport {
    let days = num_days.clamp(1, 60);
    let mut entries = Vec::with_capacity(days);

    // Floor to local 0h UT + 0.5 day step (`moonrise.c` iterates day by day)
    let jd0 = (start_jd_ut - 0.5).floor() + 0.5;

    for i in 0..days {
        let jd_ut = jd0 + (i as f64);
        let (y, m, d, _, _, _, _) = jtocal_parts(jd_ut);
        let y_dec = (y as f64) + ((m as f64) - 1.0) / 12.0 + (d as f64) / 365.25;
        let dt_sec = deltat(y_dec);
        let jd_tdt = jd_ut + dt_sec / 86400.0;

        let dm = domoon(
            jd_tdt,
            jd_ut,
            obs_lon_deg,
            obs_lat_deg,
            obs_height_m,
            temp_c,
            pressure_mb,
        );

        entries.push(MoonriseDayEntry {
            day_index: i + 1,
            date_str: format!("{:04}-{:02}-{:02}", y, m, d),
            jd_ut_noon: jd_ut,
            rise_ut: if dm.rise_transit_set.has_rise {
                Some(dm.rise_transit_set.rise_hms_ut.clone())
            } else {
                None
            },
            transit_ut: if dm.rise_transit_set.has_transit {
                Some(dm.rise_transit_set.transit_hms_ut.clone())
            } else {
                None
            },
            set_ut: if dm.rise_transit_set.has_set {
                Some(dm.rise_transit_set.set_hms_ut.clone())
            } else {
                None
            },
            rise_jd: dm.rise_transit_set.jd_rise_ut,
            transit_jd: dm.rise_transit_set.jd_transit_ut,
            set_jd: dm.rise_transit_set.jd_set_ut,
            illuminated_pct: dm.illuminated_fraction * 100.0,
            phase_description: dm.phase_description,
            app_ra_hms: dm.app_ra_hms,
            app_dec_dms: dm.app_dec_dms,
            horiz_parallax_dms: dm.horiz_parallax_dms,
        });
    }

    MoonriseReport {
        start_jd_ut: jd0,
        num_days: days,
        obs_lon_deg,
        obs_lat_deg,
        obs_height_m,
        entries,
    }
}

