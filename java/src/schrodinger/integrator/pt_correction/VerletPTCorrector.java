package schrodinger.integrator.pt_correction;

import schrodinger.QuantumLevel;

/**
 * Perturbative correction for the Verlet integration method.
 */
public class VerletPTCorrector implements PTCorrector {

    @Override
    public void computeAndSet(QuantumLevel level) {
        double result = 0;
        for (int i = 0; i < level.getGrid().getNumberOfPoints(); i++) {
            result += Math.pow(level.QTildeAtGridPoint(i) * level.psi[i], 2);
        }
        double h = level.getGrid().getStepSizeYCoordinate();
        result = result * h * h * h / (24.0 * level.getMass());
        level.perturbativeCorrectionToEnergy = result;
    }
}