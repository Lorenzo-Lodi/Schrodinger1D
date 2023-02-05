package net.tinvention.schrodinger.integrator;

import net.tinvention.schrodinger.potential.DressedPotential;

public interface Integrator {
	public double propagate(double x, double y0, double y1, double stepSize, DressedPotential potential);
}
