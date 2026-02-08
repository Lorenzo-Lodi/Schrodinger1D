package schrodinger.integrator;

import schrodinger.EnergyLevel;
import schrodinger.potential.SchrodingerSystem;

public interface Integrator {
    double propagateForward(double[] psi, int n, SchrodingerSystem potential);

    double propagateBackward(double[] psi, int n, SchrodingerSystem potential);

    // TODO computePerturbativeCorrection should probably go into a separate interface. Consider refactoring.
    default double computePerturbativeCorrection(EnergyLevel level, SchrodingerSystem qTilde) {
        return 0;
    }

}
