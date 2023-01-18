package net.tinvention.schrodinger.integrator;

import net.tinvention.schrodinger.potential.DressedPotential;

public class TaylorThreePoints implements Integrator {
	public double propagate(double x, double y0, double y1, double h, DressedPotential v) {
		return y1 * (2.0d + h * h * v.value(x)) - y0;
	}
}
