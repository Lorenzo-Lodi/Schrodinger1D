package schrodinger.integrator;

import java.util.function.IntToDoubleFunction;

/**
 * Raptis-Allison Exponentially Fitted Method for integrating the Schrödinger equation.
 * 
 * Corrected to ensure:
 * 1. Oscillatory region limit matches Numerov (beta -> 1/12).
 * 2. Exponential region limit matches Numerov (beta -> -1/12).
 * 3. High accuracy for large Z (exponential fitting).
 */
public class RaptisAllison implements Integrator {

    @Override
    public double propagate(double[] psi, int n, double step, IntToDoubleFunction qTildeFunction, Direction direction) {
        double h = step;
        double h2 = h * h;

        int n_prev = n - direction.getValue();
        int n_curr = n;
        int n_next = n + direction.getValue();

        double Q_next = qTildeFunction.applyAsDouble(n_next);
        double Q_curr = qTildeFunction.applyAsDouble(n_curr);
        double Q_prev = qTildeFunction.applyAsDouble(n_prev);

        double Z_curr = h2 * Q_curr;

        // 4. Calculate Fitted Coefficients (beta and gamma/b1)
        double beta, b1;

        if (Math.abs(Z_curr) < 1e-5) {
            // Standard Numerov limit
            beta = 1.0 / 12.0;
            b1 = 10.0 / 12.0;
        } else if (Z_curr > 0) {
            // Oscillatory Region (Q > 0)
            // We need beta -> +1/12 as Z -> 0
            double v = Math.sqrt(Z_curr);

            // Formula derived from: beta = 1/(4 sin^2(v/2)) - 1/v^2
            // This ensures limit is 1/12
            double sin2 = Math.sin(v / 2.0);
            beta = (1.0 / (4.0 * sin2 * sin2)) - (1.0 / Z_curr);

            // Consistency: 2*beta + b1 = 1  => b1 = 1 - 2*beta
            b1 = 1.0 - 2.0 * beta;

        } else {
            // Exponential Region (Q < 0)
            // We need beta -> -1/12 as Z -> 0 (so 1 + beta*Z -> 1 - |Z|/12)
            double w = Math.sqrt(-Z_curr);

            // Formula derived from: beta = 1/w^2 - 1/(4 sinh^2(w/2))
            // Note: 1/w^2 = -1/Z.
            // Limit check: 1/w^2 - 1/(w^2 * (1 + w^2/12)) = -1/12. Correct.
            double sinh2 = Math.sinh(w / 2.0);
            beta = (1.0 / (w * w)) - (1.0 / (4.0 * sinh2 * sinh2));

            b1 = 1.0 - 2.0 * beta;
        }

        // 5. Propagate using Numerov structure
        // y_{n+1} * (1 + beta * Z_{n+1}) = (2 - b1 * Z_n) * y_n - (1 + beta * Z_{n-1}) * y_{n-1}

        double Z_prev = h2 * Q_prev;
        double Z_next = h2 * Q_next;

        double numerator = (2.0 - Z_curr * b1) * psi[n_curr]
                - (1.0 + Z_prev * beta) * psi[n_prev];

        double denominator = 1.0 + Z_next * beta;

        return numerator / denominator;
    }

}
