//// ---- TRAITS -------------------------
pub trait Update{
	fn update(&mut self);
	fn get_time(&self)->f64;
}

pub trait Plot{
	fn plot(&self);
}

pub trait GetDeltaT{
	fn get_dt(&self) -> f64;
}

/// ---- Functions ---------------------
pub fn get_dt_neumann(dx: f64, safety: f64) -> f64 {
	let k: f64 = 1.0;
	return safety*(dx*dx)/(2.0*k)
}

pub fn iterate_to_max_time<S>(sim: &mut S,max_time: f64)
where
	S: Update
{
	loop{
			sim.update();
			println!("Time: {:.5}", sim.get_time());
			if sim.get_time()>max_time { break; }
		}
}