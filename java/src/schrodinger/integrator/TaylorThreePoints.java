package schrodinger.integrator;

import schrodinger.QuantumState;
import schrodinger.potential.SchrodingerSystem;

public class TaylorThreePoints implements Integrator {

    public double propagate(double[] psi, int n, QuantumState state, Direction direction) {
        double hy = state.getGrid().getStepSizeYCoordinate();
        return psi[n] * (2.0d - hy * hy * state.QTildeValueAt(n)) - psi[n - direction.getValue()];
    }

    @Override
    public void computePerturbativeCorrection(QuantumState level) {
        double result = 0;
        for (int i = 0; i < level.getGrid().getNumberOfPoints(); i++) {
            result += Math.pow(level.QTildeValueAt(i) * level.psi[i], 2);
        }
        double h = level.getGrid().getStepSizeYCoordinate();
        result = result * h * h * h / (24.0 * level.getMass());
        level.perturbativeCorrectionToEnergy = result;
    }

}
