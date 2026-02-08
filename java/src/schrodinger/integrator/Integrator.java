package schrodinger.integrator;

import schrodinger.EnergyLevel;
import schrodinger.potential.TransformedQFunction;

public interface Integrator {
	double propagateForward(double[] psi, int n, TransformedQFunction potential);

	default double computePerturbativeCorrection(EnergyLevel level, TransformedQFunction qTilde) {
		return 0;
	}

}
