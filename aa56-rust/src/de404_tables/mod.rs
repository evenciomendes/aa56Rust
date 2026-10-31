//! DE404 Ephemeris Tables (`mer404`..`plu404`, `ear404`, `moonlr`, `moonlat`)
//! Converted 100% from Stephen L. Moshier's `aa-56` (`*404.c` & `plantbl.h`).

pub struct PlanTblF64 {
    pub maxargs: i8,
    pub max_harmonic: [i8; 18],
    pub max_power_of_t: i8,
    pub arg_tbl: &'static [i8],
    pub lon_tbl: &'static [f64],
    pub lat_tbl: &'static [f64],
    pub rad_tbl: &'static [f64],
    pub distance: f64,
    pub timescale: f64,
    pub trunclvl: f64,
}

pub struct PlanTblI32 {
    pub maxargs: i8,
    pub max_harmonic: [i8; 18],
    pub max_power_of_t: i8,
    pub arg_tbl: &'static [i8],
    pub lon_tbl: &'static [i32],
    pub lat_tbl: &'static [i32],
    pub rad_tbl: &'static [i32],
    pub distance: f64,
    pub timescale: f64,
    pub trunclvl: f64,
}

pub mod mer404;
pub mod ven404;
pub mod ear404;
pub mod mar404;
pub mod jup404;
pub mod sat404;
pub mod ura404;
pub mod nep404;
pub mod plu404;
pub mod mlr404;
pub mod mlat404;

pub use mer404::MER404;
pub use ven404::VEN404;
pub use ear404::EAR404;
pub use mar404::MAR404;
pub use jup404::JUP404;
pub use sat404::SAT404;
pub use ura404::URA404;
pub use nep404::NEP404;
pub use plu404::PLU404;
pub use mlr404::MOONLR;
pub use mlat404::MOONLAT;

