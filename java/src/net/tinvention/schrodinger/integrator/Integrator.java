package net.tinvention.schrodinger.integrator;

import net.tinvention.schrodinger.potential.DressedPotential;

public interface Integrator {
	public double propagate(double x, double psiOfNMinusOne, double psiOfN, double stepSize, DressedPotential potential);
}
