package schrodinger.integrator;

import schrodinger.QuantumState;

public class Vignoli implements Integrator {

    @Override
    public double propagate(double[] psi, int n, QuantumState state, Direction direction) {
        double h = state.getGrid().getStepSizeYCoordinate();
        double h2 = h * h;

        // 1. Identify neighbors
        int prev = n - direction.getValue();
        int curr = n;
        int next = n + direction.getValue();

        // 2. Calculate Z = h^2 * Q for each point
        // Note: Your Q is defined such that f'' = -Qf.
        // Positive Q -> Oscillatory. Negative Q -> Exponential.
        double Z_prev = h2 * state.QTildeValueAt(prev);
        double Z_curr = h2 * state.QTildeValueAt(curr);
        double Z_next = h2 * state.QTildeValueAt(next);

        // 3. Calculate Vignoli coefficients dynamically
        // We fix beta = 1/12 to maintain the Numerov structure, but tune gamma.
        double beta_prev = getBeta(Z_prev);
        double beta_next = getBeta(Z_next);
        double gamma_curr = getGamma(Z_curr, 1.0 / 12.0);

        // 4. Propagate
        // Standard Numerov: (1 + 1/12*Z_next)*y_next = (2 - 10/12*Z_curr)*y_curr - (1 + 1/12*Z_prev)*y_prev
        // Vignoli:          (1 + beta*Z_next)*y_next = (2 - gamma*Z_curr)*y_curr - (1 + beta*Z_prev)*y_prev

        double numerator = (2.0 - Z_curr * gamma_curr) * psi[curr]
                - (1.0 + Z_prev * beta_prev) * psi[prev];

        double denominator = 1.0 + Z_next * beta_next;

        return numerator / denominator;
    }

    /**
     * For standard "Exponentially Fitted Numerov", beta is often kept at 1/12
     * while gamma is tuned to satisfy the frequency condition.
     */
    private double getBeta(double Z) {
        return 1.0 / 12.0;
    }

    /**
     * Calculates gamma such that the method is exact for the local frequency.
     *
     * @param Z    The local value of h^2 * Q
     * @param beta The beta coefficient (usually 1/12)
     */
    private double getGamma(double Z, double beta) {
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


