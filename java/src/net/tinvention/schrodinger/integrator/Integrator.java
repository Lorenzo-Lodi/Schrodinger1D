package net.tinvention.schrodinger.integrator;

import net.tinvention.schrodinger.potential.Potential;

public interface Integrator {
	public double propagate(double x, double y0, double y1, double step, double mass, Potential v, double e);
}
