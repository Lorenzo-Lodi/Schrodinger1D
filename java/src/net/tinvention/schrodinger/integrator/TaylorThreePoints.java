package net.tinvention.schrodinger.integrator;

import net.tinvention.schrodinger.potential.Potential;

public class TaylorThreePoints implements Integrator {

	public double propagate(double x, double y0, double y1, double step, double mass, Potential v, double e) {
		return y1 * (2.0d + step * step * v.dressedValue(x, mass, e)) - y0;
	}
}
