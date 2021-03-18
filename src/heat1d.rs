extern crate dmsuite;
extern crate ndarray;
extern crate ndarray_linalg;
use dmsuite::*;
use ndarray::*;
use ndarray_linalg::*;
use crate::utils::*;

//// ---- Heat1d -------------------------
pub struct Heat1D {
    delta_t: f64,
    time: f64,
    pub u: Array1<f64>,
    pub x: Array1<f64>,
    rhs: Array2<f64>,
}

impl Heat1D {
	pub fn new(dim: usize, method: String) -> Heat1D {
		// Check Method
		assert!( method.contains("Implicit") ||  
			method.contains("Explicit") );
		// Set solution vector
		let u = Array1::<f64>::zeros( dim ); 
		// Get Fourier differentiation matrix
		let (x,d2) = fourdif(dim,2); 
		// Get timestep size
		let safety = if method.contains("Implicit") {
			20.0 } else { 0.4 };
		let dt = get_dt_neumann(x[1]-x[0],safety); 
		// Set rhs
		let m = if method.contains("Implicit") {
			let mut a =  Array2::<f64>::eye(dim)-dt*d2;
			let eye = Array2::<f64>::eye(dim);
			a.slice_mut(s![0,..])
				.assign(&eye.slice(s![0,..])); // BC LEFT
			a.slice_mut(s![dim-1,..])
				.assign(&eye.slice(s![dim-1,..])); // BC RIGHT
			a.inv().unwrap()
		} else {
			Array2::<f64>::eye(dim)+dt*d2
		};
		Heat1D { delta_t: dt, time: 0.0, u: u, x: x, rhs:m}
	}

	pub fn iterate(&mut self, max_time:f64) {
		iterate_to_max_time(self,max_time);
	}

}

impl Update for Heat1D {
	fn update(&mut self) {
			let n = self.u.len();
			// unew = (1+k*dt*D2)*uold
			self.u = self.rhs.dot(&self.u);
			// Set edge
			self.u[0] = 0.0;
			self.u[n-1] = 0.0;
			// Update time
			self.time += self.delta_t;
		}
	fn get_time(&self) -> f64{
		self.time
	}
}

use gnuplot::{AutoOption::Fix, AxesCommon, Color, Figure};
impl Plot for Heat1D {
	fn plot(&self) {
		let x: Vec<f64> = self.x.to_vec();
		let y: Vec<f64> = self.u.to_vec();
		let mut fg = Figure::new();
		let axes = fg.axes2d();
		// Scale the x data so it lines up with the width given.
		axes.lines(&x, &y, &[Color("red")]);
		// Set the window.
		axes.set_x_range(Fix(0.), Fix(x[x.len()-1]));
		axes.set_y_range(Fix(0.), Fix(1.));
		// "wxt" is the GUI name and "768,768" is the size of it in pixels.
		fg.set_terminal("wxt size 768,768", "");
		fg.show().unwrap();
	}
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_heat1d_explicit (){
        // Initialize
        let sim_n = 80;
    	let sim_time = 1.0;
		let mut sim = Heat1D::new(sim_n,String::from("Explicit"));
		// Initial disturbance
		sim.u = sim.x.mapv(|x| (0.5*x).sin()); 
		sim.iterate(sim_time);
		// sim.plot();
        // Check
        let a = sim.u[sim_n/2];
        assert!(  0.75 < a && 0.8 > a );
    }

    #[test]
    fn test_heat1d_implicit (){
        // Initialize
        let sim_n = 80;
    	let sim_time = 1.0;
		let mut sim = Heat1D::new(sim_n,String::from("Implicit"));
		// Initial disturbance
		sim.u = sim.x.mapv(|x| (0.5*x).sin()); 
		sim.iterate(sim_time);
		// sim.plot();
        // Check
        let a = sim.u[sim_n/2];
        assert!(  0.75 < a && 0.8 > a );
    }

}