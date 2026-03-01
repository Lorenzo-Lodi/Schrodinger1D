package schrodinger.pt_correction;

import schrodinger.QuantumState;

/**
 * Perturbative correction for the Taylor Three Points integration method.
 */
public class TaylorThreePointsPTCorrection implements PTCorrection {

    @Override
    public void compute(QuantumState level) {
        double result = 0;
        for (int i = 0; i < level.getGrid().getNumberOfPoints(); i++) {
            result += Math.pow(level.QTildeAtGridPoint(i) * level.psi[i], 2);
        }
        double h = level.getGrid().getStepSizeYCoordinate();
        result = result * h * h * h / (24.0 * level.getMass());
        level.perturbativeCorrectionToEnergy = result;
    }
}