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
        level.perturbativeCorrectionToEnergy = 0;
    }

}
