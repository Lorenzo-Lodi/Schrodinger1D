package net.tinvention.schrodinger.integrator;

import net.tinvention.schrodinger.potential.DressedPotential;

public class TaylorThreePoints implements Integrator {
	public double propagate(double yOfN, double psiOfNMinusOne, double psiOfN, double stepSize, DressedPotential v) {
		return psiOfN * (2.0d - stepSize * stepSize * v.value(yOfN)) - psiOfNMinusOne;
	}
}
