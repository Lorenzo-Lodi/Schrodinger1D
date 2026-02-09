package schrodinger.integrator;

import schrodinger.EnergyLevel;
import schrodinger.potential.SchrodingerSystem;

public class TaylorThreePoints implements Integrator {

    public double propagate(double[] psi, int n, SchrodingerSystem system, Direction direction) {
        double y = system.getGrid().getYValue(n);
        double stepSizeSquared = Math.pow(system.getGrid().getStepSizeYCoordinate(), 2);
        return psi[n] * (2.0d - stepSizeSquared * system.QTilde(y)) - psi[n - direction.getValue()];
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
