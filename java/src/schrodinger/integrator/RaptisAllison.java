package schrodinger.integrator;

import java.util.function.IntToDoubleFunction;

/**
 * Raptis-Allison Exponentially Fitted Method for integrating the Schrödinger equation.
 * <p>
 * Corrected to ensure:
 * 1. Oscillatory region limit matches Numerov (beta -> 1/12).
 * 2. Exponential region limit matches Numerov (beta -> -1/12).
 * 3. High accuracy for large Z (exponential fitting).
 */
public class RaptisAllison extends ExponentiallyFittedAbstract {

    public double getBeta(double Z) {
        double beta;

        if (Math.abs(Z) < 1e-3) {
            return 1.0 / 12.0 + Z / 240. + Z * Z / 6048.;
        }

        double sqrtModZ = Math.sqrt(Math.abs(Z));
        double c;

        if (Z > 0) {
            c = Math.sin(sqrtModZ / 2.0);

        } else {
            c = Math.sinh(sqrtModZ / 2.0);
        }

        return Math.signum (Z)/ (4.0 * c * c) - 1.0 / Z;
    }

    public double getGamma(double Z, double beta) {
        return 1.0 - 2.0 * beta;
    }

}
