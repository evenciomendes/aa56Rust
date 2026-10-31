//! `aa56-rust` — Complete 100% Rust port of Stephen L. Moshier's `aa-56` Astronomical Almanac.
//! Includes all 11 DE404 ephemeris tables (`mer404`..`plu404`, `ear404`, `mlr404`, `mlat404`),
//! `gplan`, `kepler`, `domoon`, `conjunct`, `moonrise`, `sun`, `rplanet`, `rstar`, `fk4fk5`,
//! `trnsit`, `altaz`, `constel`, `deltat`, and embedded astronomical catalogs (`star.cat`,
//! `star1950.cat`, `messier.cat`, `orbit.cat`, `aa.ini`).

pub mod conjunct;
pub mod constel;
pub mod de404_tables;
pub mod domoon;
pub mod earth_orientation;
pub mod gplan;
pub mod kepler;
pub mod math_time;
pub mod moonrise;
pub mod reductions;
pub mod stars_catalog;
pub mod transit;

/// Embedded `aa.ini` configuration from `aa-56`.
pub const CATALOG_AA_INI: &str = include_str!("../catalogs/aa.ini");
/// Embedded `star.cat` FK5/FK4 star catalog from `aa-56`.
pub const CATALOG_STAR: &str = include_str!("../catalogs/star.cat");
/// Embedded `star1950.cat` FK4 B1950 star catalog.
pub const CATALOG_STAR_1950: &str = include_str!("../catalogs/star1950.cat");
/// Embedded `messier.cat` Messier deep-sky catalog from `aa-56`.
pub const CATALOG_MESSIER: &str = include_str!("../catalogs/messier.cat");
/// Embedded `orbit.cat` comets and asteroids catalog from `aa-56`.
pub const CATALOG_ORBIT: &str = include_str!("../catalogs/orbit.cat");

