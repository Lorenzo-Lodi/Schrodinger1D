package schrodinger.integrator;

import schrodinger.QuantumState;
import schrodinger.potential.SchrodingerSystem;

public class TaylorThreePoints implements Integrator {

    public double propagate(double[] psi, int n, SchrodingerSystem system, Direction direction) {
        double stepSizeSquared = Math.pow(system.getGrid().getStepSizeYCoordinate(), 2);
        return psi[n] * (2.0d - stepSizeSquared * system.QTildeValueAt(n)) - psi[n - direction.getValue()];
    }

    @Override
    public double computePerturbativeCorrection(QuantumState level, SchrodingerSystem system) {
        double result = 0;
        if (!level.isPsiNormalized()) {
            level.normalizePsi();
        }
        for (int i = 0; i < system.getGrid().getNumberOfPoints(); i++) {
            result += Math.pow(system.QTildeValueAt(i) * level.psi[i], 2);
        }
        double h = system.getGrid().getStepSizeYCoordinate();
        result = result * h * h * h / (24.0 * system.getMass());
        return result;
    }

}
