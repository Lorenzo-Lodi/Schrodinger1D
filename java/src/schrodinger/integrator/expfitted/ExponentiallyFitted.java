package schrodinger.integrator.expfitted;

public class ExponentiallyFitted extends ExponentiallyFittedAbstract {

    /**
     * Computes b0 (beta) for S1 EF Numerov.
     * Exact formula: b0 = [Z - 2(1 - c)] / [2 Z (1 - c)], where c = cos(sqrt|Z|) or cosh(sqrt|Z|)
     */
    public double getBeta(double Z) {

        if (Math.abs(Z) < 0.01) {
            double Z2 = Z * Z;
            double Z3 = Z * Z2;
            return 1.0 / 12.0 + Z / 240. + Z2 / 6048. + Z3 / 172800.;
        }

        double sqrtModZ = Math.sqrt(Math.abs(Z));
        double c;

        if (Z > 0) {
            c = Math.sin(sqrtModZ / 2.0);

        } else {
            c = Math.sinh(sqrtModZ / 2.0);
        }

        return Math.signum(Z) / (4.0 * c * c) - 1.0 / Z;

// The following code is mathematically equivalent, but leads to slightly more numerical error for small Z
//  (although we use Taylor expansion for small Z anyway).
//        if (Z > 0) {  // Oscillatory Region (Q > 0)
//            c = 1.0 - Math.cos(sqrtModZ);
//        } else { // Exponential Region (Q < 0)
//            c = 1.0 - Math.cosh(sqrtModZ);
//        }
//
//        return (Z - 2.0 * c) / (2.0 * Z * c);
    }

    public double getGamma(double Z, double beta) {
        return 1.0 - 2.0 * beta;
    }
}
