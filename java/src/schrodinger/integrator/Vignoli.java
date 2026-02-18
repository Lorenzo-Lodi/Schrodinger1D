package schrodinger.integrator;

public class Vignoli extends ExponentiallyFittedAbstract {

    double getBeta(double Z) {
        return 1.0 / 12.0;
    }

    double getGamma(double Z, double beta) {

        // Limit for small Z to avoid division by zero
        if (Math.abs(Z) < 1e-4) {
            return 5.0 / 6.0; // Standard Numerov value
        }

        if (Z > 0) { // Oscillatory Region (Q > 0)
            double sqrtZ = Math.sqrt(Z);
            double cosSqrtZ = Math.cos(sqrtZ);
            return (2.0 - 2.0 * cosSqrtZ * (1.0 + beta * Z)) / Z;
        } else { // Exponential Region (Q < 0)
            double sqrtModZ = Math.sqrt(-Z);
            double coshSqrtZ = Math.cosh(sqrtModZ);
            return (2.0 - 2.0 * coshSqrtZ * (1.0 + beta * Z)) / Z;
        }
    }
}


