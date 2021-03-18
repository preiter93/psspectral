
use crate::utils::Plot;
use crate::heat1d::Heat1D;

mod heat1d;
mod utils;

fn main() {
    test_heat1d()
}

// Amount of grid points we use
const SIM_N: usize = 180;

// Amount of time we will simulate
const SIM_TIME: f64 = 10.;

pub fn test_heat1d() {
	// Initialize
	let mut sim = Heat1D::new(SIM_N,String::from("Implicit"));
	sim.u = sim.x.mapv(|x| (0.5*x).sin()); // Initial disturbance
	sim.iterate(SIM_TIME);
	sim.plot();
	println!("{:?}", sim.u[SIM_N/2]);
}