package schrodinger.integrator;

import schrodinger.EnergyLevel;
import schrodinger.potential.SchrodingerSystem;

public class TaylorThreePoints implements Integrator {

    @Override
    public double propagateForward(double[] psi, int n, SchrodingerSystem qTilde) {
        return propagate(psi, n, qTilde, 1);
    }

    @Override
    public double propagateBackward(double[] psi, int n, SchrodingerSystem qTilde) {
        return propagate(psi, n, qTilde, -1);
    }

    private double propagate(double[] psi, int n, SchrodingerSystem qTilde, int direction) {
        double y = qTilde.getGrid().getYValue(n);
        double stepSizeSquared = Math.pow(qTilde.getGrid().getStepSizeYCoordinate(), 2);
        return psi[n] * (2.0d - stepSizeSquared * qTilde.QTilde(y)) - psi[n - direction];
    }

    @Override
    public double computePerturbativeCorrection(EnergyLevel level, SchrodingerSystem qTilde) {
        double result = 0;
        for (int i = 0; i < qTilde.getGrid().getNumberOfPoints(); i++) {
            result += Math.pow(qTilde.QTilde(qTilde.getGrid().getYValue(i)) * level.psi[i], 2);
        }
        result = result * Math.pow(qTilde.getGrid().getStepSizeYCoordinate(), 3) / (24 * qTilde.getMass());
        return result;
    }

}
