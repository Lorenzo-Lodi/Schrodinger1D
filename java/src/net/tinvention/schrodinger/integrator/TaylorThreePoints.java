package net.tinvention.schrodinger.integrator;

import net.tinvention.schrodinger.EnergyLevel;
import net.tinvention.schrodinger.potential.DressedPotential;

public class TaylorThreePoints implements Integrator {

	@Override
	public double propagate(double yOfN, double psiOfNMinusOne, double psiOfN, int n, double stepSize,
			DressedPotential qTilde) {
		return psiOfN * (2.0d - stepSize * stepSize * qTilde.value(yOfN)) - psiOfNMinusOne;
	}

	@Override
	public double computePerturbativeCorrection(EnergyLevel level, DressedPotential qTilde) {
		double result = 0;
		for (int i = 0; i < qTilde.getGrid().getNumberOfPoints(); i++) {
			result += Math.pow(qTilde.value(qTilde.getGrid().getYValue(i)) * level.psi[i], 2);
		}
		result = result * Math.pow(qTilde.getGrid().getStepSizeYCoordinate(), 2) / (24 * qTilde.getMass());
		return result;
	}

}
