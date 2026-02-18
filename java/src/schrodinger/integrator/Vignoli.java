package schrodinger.integrator;

public class Vignoli extends ExponentiallyFittedAbstract {

    /**
     * For standard "Exponentially Fitted Numerov", beta is often kept at 1/12
     * while gamma is tuned to satisfy the frequency condition.
     */
    double getBeta(double Z) {
        return 1.0 / 12.0;
    }

    /**
     * Calculates gamma such that the method is exact for the local frequency.
     *
     * @param Z    The local value of h^2 * Q
     * @param beta The beta coefficient (usually 1/12)
     */
    double getGamma(double Z, double beta) {
        // Limit for small Z to avoid division by zero
        if (Math.abs(Z) < 1e-4) {
            return 5.0 / 6.0; // Standard Numerov value
        }

        if (Z > 0) {
            // Oscillatory Region (Q > 0)
            // Condition: y_{n+1} + y_{n-1} = 2*cos(sqrt(Z)) * y_n
            // Plugging into recurrence yields:
            double sqrtZ = Math.sqrt(Z);
            double cosSqrtZ = Math.cos(sqrtZ);

            // Formula derived from: 2*cos = (2 - gamma*Z) / (1 + beta*Z)
            // Rearranging for gamma:
            return (2.0 - 2.0 * cosSqrtZ * (1.0 + beta * Z)) / Z;

        } else {
            // Exponential Region (Q < 0)
            // Condition: y_{n+1} + y_{n-1} = 2*cosh(sqrt(-Z)) * y_n
            // We use -Z because Z is negative here.
            double sqrtModZ = Math.sqrt(-Z);
            double coshSqrtZ = Math.cosh(sqrtModZ);

            // Formula derived from: 2*cosh = (2 - gamma*Z) / (1 + beta*Z)
            // Rearranging for gamma:
            return (2.0 - 2.0 * coshSqrtZ * (1.0 + beta * Z)) / Z;
        }
    }
}


