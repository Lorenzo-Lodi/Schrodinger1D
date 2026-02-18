package schrodinger.integrator;

public class ExponentiallyFitted extends ExponentiallyFittedAbstract {

    /**
     * Computes b0 (beta) for S1 EF Numerov.
     * Exact formula: b0 = [Z - 2(1 - c)] / [2 Z (1 - c)], where c = cos(sqrt|Z|) or cosh(sqrt|Z|)
     */
    double getBeta(double Z) {
        if (Math.abs(Z) < 1e-10) {
            return 1.0 / 12.0 + Z / 240.;  // Standard Numerov value
        }

        double sqrtZ = Math.sqrt(Math.abs(Z));
        double denom;

        if (Z > 0) {  // Oscillatory Region (Q > 0)
            denom = 1.0 - Math.cos(sqrtZ);
        } else { // Exponential Region (Q < 0)
            denom = 1.0 - Math.cosh(sqrtZ);
        }

        return (Z - 2.0 * denom) / (2.0 * Z * denom);
    }

    double getGamma(double Z, double beta) {
        return 1.0 - 2.0 * beta;
    }
}
