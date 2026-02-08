package schrodinger.integrator;

import schrodinger.EnergyLevel;
import schrodinger.potential.TransformedQFunction;

public interface Integrator {
    double propagateForward(double[] psi, int n, TransformedQFunction potential);

    double propagateBackward(double[] psi, int n, TransformedQFunction potential);

    // TODO computePerturbativeCorrection should probably go into a separate interface. Consider refactoring.
    default double computePerturbativeCorrection(EnergyLevel level, TransformedQFunction qTilde) {
        return 0;
    }

}
