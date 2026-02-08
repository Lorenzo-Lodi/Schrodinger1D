package schrodinger.integrator;

import schrodinger.EnergyLevel;
import schrodinger.potential.SchrodingerSystem;

public class TaylorThreePoints implements Integrator {

    @Override
    public double propagateForward(double[] psi, int n, SchrodingerSystem system) {
        return propagate(psi, n, system, 1);
    }

    @Override
    public double propagateBackward(double[] psi, int n, SchrodingerSystem system) {
        return propagate(psi, n, system, -1);
    }

    private double propagate(double[] psi, int n, SchrodingerSystem system, int direction) {
        double y = system.getGrid().getYValue(n);
        double stepSizeSquared = Math.pow(system.getGrid().getStepSizeYCoordinate(), 2);
        return psi[n] * (2.0d - stepSizeSquared * system.QTilde(y)) - psi[n - direction];
    }

    @Override
    public double computePerturbativeCorrection(EnergyLevel level, SchrodingerSystem system) {
        double result = 0;
        for (int i = 0; i < system.getGrid().getNumberOfPoints(); i++) {
            result += Math.pow(system.QTilde(system.getGrid().getYValue(i)) * level.psi[i], 2);
        }
        result = result * Math.pow(system.getGrid().getStepSizeYCoordinate(), 3) / (24 * system.getMass());
        return result;
    }

}
