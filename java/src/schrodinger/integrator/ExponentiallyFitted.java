package schrodinger.integrator;

public class ExponentiallyFitted extends ExponentiallyFittedAbstract {

    /**
     * Computes b0 (beta) for S1 EF Numerov.
     * Exact formula: b0 = [Z - 2(1 - c)] / [2 Z (1 - c)], where c = cos(sqrt|Z|) or cosh(sqrt|Z|)
     */
    public double getBeta(double Z) {

        if (Math.abs(Z) < 1e-3) {
            return 1.0 / 12.0 + Z / 240. + Z * Z / 6048.;  // Standard Numerov value
        }

        double sqrtModZ = Math.sqrt(Math.abs(Z));
        double c;

        if (Z > 0) {  // Oscillatory Region (Q > 0)
            c = 1.0 - Math.cos(sqrtModZ);
        } else { // Exponential Region (Q < 0)
            c = 1.0 - Math.cosh(sqrtModZ);
        }

        return (Z - 2.0 * c) / (2.0 * Z * c);
    }

    public double getGamma(double Z, double beta) {
        return 1.0 - 2.0 * beta;
    }
}
