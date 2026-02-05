package net.tinvention.schrodinger.integrator;

import net.tinvention.schrodinger.EnergyLevel;
import net.tinvention.schrodinger.potential.DressedPotential;

public interface Integrator {
	double propagateForward(double[] psi, int n, DressedPotential potential);

	default double computePerturbativeCorrection(EnergyLevel level, DressedPotential qTilde) {
		return 0;
	}

}
