package schrodinger.integrator;

import schrodinger.EnergyLevel;
import schrodinger.potential.SchrodingerSystem;

public interface Integrator {
    double propagateForward(double[] psi, int n, SchrodingerSystem system);

    double propagateBackward(double[] psi, int n, SchrodingerSystem system);

    // TODO computePerturbativeCorrection should probably go into a separate interface. Consider refactoring.
    default double computePerturbativeCorrection(EnergyLevel level, SchrodingerSystem system) {
        return 0;
    }

}
