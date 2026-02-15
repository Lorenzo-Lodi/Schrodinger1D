package schrodinger.integrator;

import schrodinger.QuantumState;

public class Numerov implements Integrator {

    public double propagate(double[] psi, int n, QuantumState state, Direction direction) {
        double hy = state.getGrid().getStepSizeYCoordinate();
        return (psi[n] * (2.0d - 5. * hy * hy * state.QTildeValueAt(n) / 6.) -
                (1. + hy * hy * state.QTildeValueAt(n - direction.getValue()) / 12.)
                        * psi[n - direction.getValue()]) /
                (1. + hy * hy * state.QTildeValueAt(n + direction.getValue()) / 12.);
    }

    @Override
    public void computePerturbativeCorrection(QuantumState level) {
        double result = 0;
        for (int i = 1; i < level.getGrid().getNumberOfPoints() - 1; i++) {
            result += Math.pow(level.QTildeValueAt(i + 1) * level.psi[i + 1] - level.QTildeValueAt(i - 1) * level.psi[i - 1], 2);
        }
        result = result / 4.;
        double hy = level.getGrid().getStepSizeYCoordinate();
        result = result * hy * hy * hy / (480.0 * level.getMass());
        level.perturbativeCorrectionToEnergy = result;
    }

}
