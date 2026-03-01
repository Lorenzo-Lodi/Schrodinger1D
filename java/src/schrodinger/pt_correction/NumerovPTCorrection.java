package schrodinger.pt_correction;

import schrodinger.QuantumState;

/**
 * Perturbative correction for the Numerov integration method.
 */
public class NumerovPTCorrection implements PTCorrection {

    @Override
    public void compute(QuantumState level) {
        double result = 0;
        for (int i = 1; i < level.getGrid().getNumberOfPoints() - 1; i++) {
            result += Math.pow(level.QTildeAtGridPoint(i + 1) * level.psi[i + 1] - level.QTildeAtGridPoint(i - 1) * level.psi[i - 1], 2);
        }
        result = result / 4.;
        double hy = level.getGrid().getStepSizeYCoordinate();
        result = result * hy * hy * hy / (480.0 * level.getMass());
        level.perturbativeCorrectionToEnergy = result;
    }
}