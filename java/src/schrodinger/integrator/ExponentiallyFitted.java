package schrodinger.integrator;

import schrodinger.QuantumState;

public class ExponentiallyFitted implements Integrator {

    public double propagate(double[] psi, int n, QuantumState state, Direction direction) {
        double h = state.getGrid().getStepSizeYCoordinate();
        double h2 = h * h;

        // 1. Identify neighbors
        int prev = n - direction.getValue();
        int curr = n;
        int next = n + direction.getValue();

        // 2. Calculate Z = h^2 * Q for each point
        double Z_prev = h2 * state.QTildeValueAt(prev);
        double Z_curr = h2 * state.QTildeValueAt(curr);
        double Z_next = h2 * state.QTildeValueAt(next);

        // 3. Calculate S1 EF coefficients dynamically
        double beta_prev = getBetaS1(Z_prev);
        double beta_next = getBetaS1(Z_next);
        double beta_curr = getBetaS1(Z_curr);  // Needed for gamma_curr
        double gamma_curr = 1.0 - 2.0 * beta_curr;

        // 4. Propagate: (1 + β Z_next) y_next = (2 - γ Z_curr) y_curr - (1 + β Z_prev) y_prev
        double numerator = (2.0 - Z_curr * gamma_curr) * psi[curr]
                - (1.0 + Z_prev * beta_prev) * psi[prev];

        double denominator = 1.0 + Z_next * beta_next;

        return numerator / denominator;
    }

    /**
     * Computes b0 (beta) for S1 EF Numerov.
     * Exact formula: b0 = [Z - 2(1 - c)] / [2 Z (1 - c)], where c = cos(sqrt|Z|) or cosh(sqrt|Z|)
     */
    private double getBetaS1(double Z) {
        if (Math.abs(Z) < 1e-12) {
            return 1.0 / 12.0;  // Classical limit
        }

        double phi;
        double c;  // cos or cosh

        if (Z > 0) {
            // Oscillatory: c = cos(phi)
            phi = Math.sqrt(Z);
            c = Math.cos(phi);
        } else {
            // Exponential: c = cosh(phi)
            phi = Math.sqrt(-Z);
            c = Math.cosh(phi);
        }

        double denom_factor = 1.0 - c;
        if (Math.abs(denom_factor) < 1e-12) {
            return 1.0 / 12.0;
        }

        return (Z - 2.0 * denom_factor) / (2.0 * Z * denom_factor);
    }
}
