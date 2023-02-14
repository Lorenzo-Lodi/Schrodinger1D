package net.tinvention.schrodinger.integrator;

import net.tinvention.schrodinger.EnergyLevel;
import net.tinvention.schrodinger.potential.DressedPotential;

public interface Integrator {
	public double propagate(double psiOfNMinusOne, double psiOfN, int n, 
			DressedPotential potential);

	public default double computePerturbativeCorrection(EnergyLevel level, DressedPotential qTilde) {
		return 0;
	}

}
