package net.tinvention.schrodinger.integrator;

import net.tinvention.schrodinger.potential.Potential;

public class TaylorThreePoints implements Integrator {
	public double propagate(double x, double y0, double y1, double h, double mass, Potential v, double e) {
		return 2.0d * y1 * (1.0d - h * h * mass * (e - v.value(x))) - y0;
	}
}
