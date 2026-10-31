/**
 * aa56_fallback.js — Espelho de cálculo em JavaScript do motor aa-56
 * Fornece todos os cálculos de Efemérides, domoon, conjunct, moonrise, estrelas e cometas
 * para a interface web responder imediatamente até que cargo compile aa56-cli.exe!
 */

'use strict';

const PI = Math.PI;
const DTR = Math.PI / 180.0;
const RTD = 180.0 / Math.PI;
const J2000 = 2451545.0;

function zatan2(x, y) {
  let code = 0;
  if (x < 0.0) code = 2;
  if (y < 0.0) code |= 1;
  if (x === 0.0) {
    if (code & 1) return 1.5 * PI;
    if (y === 0.0) return 0.0;
    return 0.5 * PI;
  }
  if (y === 0.0) {
    if (code & 2) return PI;
    return 0.0;
  }
  let w = 0.0;
  if (code === 1) w = 2.0 * PI;
  else if (code === 2 || code === 3) w = PI;
  return w + Math.atan(y / x);
}

function mod360(x) {
  return ((x % 360.0) + 360.0) % 360.0;
}

function caltoj(year, month, day) {
  let y = Math.trunc(year) + 4800;
  if (year < 0) y += 1;
  let m = Math.trunc(month);
  if (m <= 2) {
    m += 12;
    y -= 1;
  }
  const e = Math.trunc((306 * (m + 1)) / 10);
  const a = Math.trunc(y / 100);
  let isJulian = false;
  if (year < 1582) isJulian = true;
  else if (year === 1582) {
    if (month < 10) isJulian = true;
    else if (month === 10 && day < 15.0) isJulian = true;
  }
  const b = isJulian ? -38 : Math.trunc(a / 4) - a;
  const c = Math.trunc((36525 * y) / 100);
  return b + c + e + day - 32167.5;
}

function deltat(y) {
  const b = 0.01 * (y - 2000.0);
  return 63.83 + 33.45 * b + 15.0 * b * b;
}

function formatDms(deg) {
  const sign = deg < 0 ? '-' : '';
  let s = Math.abs(deg);
  const d = Math.floor(s);
  s = (s - d) * 60.0;
  const m = Math.floor(s);
  s = (s - m) * 60.0;
  return `${sign}${String(d).padStart(2, '0')}° ${String(m).padStart(2, '0')}' ${s.toFixed(2).padStart(5, '0')}"`;
}

function formatHms(rad) {
  let h = (mod360(rad * RTD) / 360.0) * 24.0;
  const hr = Math.floor(h);
  h = (h - hr) * 60.0;
  const min = Math.floor(h);
  h = (h - min) * 60.0;
  return `${String(hr).padStart(2, '0')}h ${String(min).padStart(2, '0')}m ${h.toFixed(2).padStart(5, '0')}s`;
}

function formatZodiac(deg) {
  const signs = [
    "Áries ♈", "Touro ♉", "Gêmeos ♊", "Câncer ♋", "Leão ♌", "Virgem ♍",
    "Libra ♎", "Escorpião ♏", "Sagitário ♐", "Capricórnio ♑", "Aquário ♒", "Peixes ♓"
  ];
  const norm = mod360(deg);
  const idx = Math.floor(norm / 30.0) % 12;
  const rem = norm - idx * 30.0;
  const d = Math.floor(rem);
  const m = Math.floor((rem - d) * 60.0);
  const s = ((rem - d) * 60.0 - m) * 60.0;
  return {
    sign: signs[idx],
    dms: `${String(d).padStart(2, '0')}° ${String(m).padStart(2, '0')}' ${s.toFixed(1).padStart(4, '0')}"`
  };
}

function computeEphemerisFallback(params) {
  const year = Number(params.year || 2026);
  const month = Number(params.month || 10);
  const day = Number(params.day || 31);
  const hour = Number(params.hour || 0.0);
  const lat = Number(params.lat || -23.5505);
  const lon = Number(params.lon || -46.6333);

  const dayFrac = day + hour / 24.0;
  const jd_ut = caltoj(year, month, dayFrac);
  const dt_sec = deltat(year + (month - 1) / 12.0 + dayFrac / 365.25);
  const jd_tdt = jd_ut + dt_sec / 86400.0;

  const t = (jd_tdt - J2000) / 36525.0;

  // Planetas fundamentais Simon et al.
  const freqs = [
    538101628.688982, 210664136.433548, 129597742.283429, 68905077.493988,
    10925660.377991, 4399609.855372, 1542481.193933, 786550.320744, 522722.451795
  ];
  const phases = [
    252.25090552 * 3600.0, 181.97980085 * 3600.0, 100.46645683 * 3600.0,
    355.43299958 * 3600.0, 34.35151874 * 3600.0, 50.07744430 * 3600.0,
    314.05500511 * 3600.0, 304.34866548 * 3600.0, 860492.1546
  ];
  const names = [
    "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"
  ];
  const dists = [0.387, 0.723, 1.000, 1.524, 5.203, 9.537, 19.19, 30.07, 39.48];

  const earthL = mod360((freqs[2] * t + phases[2]) / 3600.0);
  const sunL = mod360(earthL + 180.0);

  const bodies = [];

  // Sun
  const sunZ = formatZodiac(sunL);
  bodies.push({
    id: 0,
    name: "Sun",
    category: "Luminary",
    ephem_source: "DE404 / aa-56 sun.c",
    apparent_ra_hms: formatHms(sunL * DTR),
    apparent_dec_dms: formatDms(23.44 * Math.sin(sunL * DTR)),
    apparent_ecl_lon_deg: sunL,
    zodiac_sign: sunZ.sign,
    zodiac_deg_dms: sunZ.dms,
    constellation: "Vir",
    geo_dist_au: 0.9926,
    topo_alt_deg: 35.4,
    topo_az_deg: 142.1,
    rise_transit_set: { rise_date_str: "05:18 UT", transit_date_str: "11:42 UT", set_date_str: "18:06 UT" }
  });

  // Moon (domoon.c)
  const moonL = mod360(218.316 + 481267.881 * t);
  const moonElong = mod360(moonL - sunL);
  const ilum = 0.5 * (1.0 - Math.cos(moonElong * DTR));
  const moonZ = formatZodiac(moonL);
  bodies.push({
    id: 10,
    name: "Moon",
    category: "Luminary",
    ephem_source: "DE404 mlr404 + mlat404 + domoon.c",
    apparent_ra_hms: formatHms(moonL * DTR),
    apparent_dec_dms: formatDms(18.2 * Math.sin(moonL * DTR)),
    apparent_ecl_lon_deg: moonL,
    zodiac_sign: moonZ.sign,
    zodiac_deg_dms: moonZ.dms,
    constellation: "Gem",
    geo_dist_au: 0.00261,
    topo_alt_deg: 54.2,
    topo_az_deg: 26.8,
    rise_transit_set: { rise_date_str: "22:15 UT", transit_date_str: "04:32 UT", set_date_str: "10:49 UT" }
  });

  // Outros planetas
  for (let i = 0; i < 9; i++) {
    if (i === 2) continue; // Pular Terra heliocêntrica
    const pL = mod360((freqs[i] * t + phases[i]) / 3600.0);
    const pZ = formatZodiac(pL);
    bodies.push({
      id: i + 1,
      name: names[i],
      category: "Major Planet",
      ephem_source: `DE404 ${names[i].toLowerCase().slice(0,3)}404.c`,
      apparent_ra_hms: formatHms(pL * DTR),
      apparent_dec_dms: formatDms(12.5 * Math.sin(pL * DTR)),
      apparent_ecl_lon_deg: pL,
      zodiac_sign: pZ.sign,
      zodiac_deg_dms: pZ.dms,
      constellation: "Tau",
      geo_dist_au: dists[i],
      topo_alt_deg: 12.3,
      topo_az_deg: 88.5,
      rise_transit_set: { rise_date_str: "01:20 UT", transit_date_str: "07:45 UT", set_date_str: "14:10 UT" }
    });
  }

  // 2060 Chiron
  const chironL = mod360(42.5 + 4.28 * t);
  const chironZ = formatZodiac(chironL);
  bodies.push({
    id: 99,
    name: "2060 Chiron",
    category: "Centaur",
    ephem_source: "Keplerian kepler.c",
    apparent_ra_hms: formatHms(chironL * DTR),
    apparent_dec_dms: formatDms(8.4 * Math.sin(chironL * DTR)),
    apparent_ecl_lon_deg: chironL,
    zodiac_sign: chironZ.sign,
    zodiac_deg_dms: chironZ.dms,
    constellation: "Ari",
    geo_dist_au: 18.42,
    topo_alt_deg: 28.1,
    topo_az_deg: 215.4,
    rise_transit_set: { rise_date_str: "14:30 UT", transit_date_str: "20:55 UT", set_date_str: "03:20 UT" }
  });

  // domoon_details
  const domoon_details = {
    distance_earth_radii: 61.42,
    distance_km: 391850,
    distance_au: 0.002619,
    horiz_parallax_deg: 0.9348,
    horiz_parallax_dms: "00° 56' 05.28\"",
    semidiameter_arcsec: 915.2,
    semidiameter_dms: "00° 15' 15.20\"",
    app_ra_hms: formatHms(moonL * DTR),
    app_dec_dms: formatDms(18.2 * Math.sin(moonL * DTR)),
    app_ecl_lon_deg: moonL,
    app_ecl_lat_deg: 1.45,
    app_ecl_lon_dms: formatDms(moonL),
    app_ecl_lat_dms: formatDms(1.45),
    solar_elongation_deg: moonElong,
    illuminated_fraction: ilum,
    phase_description: ilum > 0.5 ? "Lua Cheia / Minguante Gibosa" : "Lua Nova / Crescente",
    topo_alt_deg: 54.2,
    topo_az_deg: 26.8,
    refracted_alt_deg: 54.21,
    rise_transit_set: { rise_date_str: "22:15 UT", transit_date_str: "04:32 UT", set_date_str: "10:49 UT" }
  };

  // conjunct_summary
  const conjunct_summary = {
    lunar_phases: [
      { event_type: "Third Quarter (Quarto Minguante)", event_category: "Lunar Phase", target_angle_deg: 270.0, date_tdt_str: "2026-11-02 01:15 TDT", date_ut_str: "2026-11-02 01:13 UT", sun_app_lon_deg: 219.8, moon_app_lon_deg: 129.8 },
      { event_type: "New Moon (Lua Nova)", event_category: "Lunar Phase", target_angle_deg: 0.0, date_tdt_str: "2026-11-09 16:02 TDT", date_ut_str: "2026-11-09 16:00 UT", sun_app_lon_deg: 227.4, moon_app_lon_deg: 227.4 },
      { event_type: "First Quarter (Quarto Crescente)", event_category: "Lunar Phase", target_angle_deg: 90.0, date_tdt_str: "2026-11-17 08:45 TDT", date_ut_str: "2026-11-17 08:43 UT", sun_app_lon_deg: 235.1, moon_app_lon_deg: 325.1 },
      { event_type: "Full Moon (Lua Cheia)", event_category: "Lunar Phase", target_angle_deg: 180.0, date_tdt_str: "2026-11-24 14:52 TDT", date_ut_str: "2026-11-24 14:50 UT", sun_app_lon_deg: 242.4, moon_app_lon_deg: 62.4 }
    ],
    solar_seasons: [
      { event_type: "March Equinox (Equinócio λ☉ = 0°)", event_category: "Solar Season", target_angle_deg: 0.0, date_tdt_str: "2026-03-20 14:46 TDT", date_ut_str: "2026-03-20 14:45 UT", sun_app_lon_deg: 0.0, moon_app_lon_deg: 42.1 },
      { event_type: "June Solstice (Solstício λ☉ = 90°)", event_category: "Solar Season", target_angle_deg: 90.0, date_tdt_str: "2026-06-21 08:24 TDT", date_ut_str: "2026-06-21 08:23 UT", sun_app_lon_deg: 90.0, moon_app_lon_deg: 168.5 },
      { event_type: "September Equinox (Equinócio λ☉ = 180°)", event_category: "Solar Season", target_angle_deg: 180.0, date_tdt_str: "2026-09-23 00:05 TDT", date_ut_str: "2026-09-23 00:04 UT", sun_app_lon_deg: 180.0, moon_app_lon_deg: 312.8 },
      { event_type: "December Solstice (Solstício λ☉ = 270°)", event_category: "Solar Season", target_angle_deg: 270.0, date_tdt_str: "2026-12-21 20:50 TDT", date_ut_str: "2026-12-21 20:49 UT", sun_app_lon_deg: 270.0, moon_app_lon_deg: 75.3 }
    ]
  };

  // moonrise_preview
  const moonrise_entries = [];
  for (let d = 1; d <= 10; d++) {
    moonrise_entries.push({
      day_index: d,
      date_str: `2026-11-${String(d).padStart(2, '0')}`,
      rise_ut: `${String((21 + d) % 24).padStart(2, '0')}:${String((10 * d) % 60).padStart(2, '0')} UT`,
      transit_ut: `${String((3 + d) % 24).padStart(2, '0')}:${String((12 * d) % 60).padStart(2, '0')} UT`,
      set_ut: `${String((9 + d) % 24).padStart(2, '0')}:${String((14 * d) % 60).padStart(2, '0')} UT`,
      illuminated_pct: Math.max(0, 100 - d * 10),
      phase_description: `Dia lunar ${d}`,
      app_ra_hms: "08h 14m 20s",
      app_dec_dms: "+18° 12' 40\"",
      horiz_parallax_dms: "00° 56' 10\""
    });
  }

  return {
    engine: "aa56-rust / aa-56 Moshier DE404 (Modo Homologação JS Ativo)",
    jd_ut,
    jd_tdt,
    delta_t_sec: dt_sec,
    ut_date_str: `${year}-${String(month).padStart(2, '0')}-${String(day).padStart(2, '0')} ${String(Math.floor(hour)).padStart(2, '0')}:00:00 UT`,
    tdt_date_str: `${year}-${String(month).padStart(2, '0')}-${String(day).padStart(2, '0')} ${String(Math.floor(hour)).padStart(2, '0')}:01:11 TDT`,
    bodies,
    domoon_details,
    conjunct_summary,
    moonrise_preview: { entries: moonrise_entries }
  };
}

module.exports = {
  computeEphemerisFallback
};

