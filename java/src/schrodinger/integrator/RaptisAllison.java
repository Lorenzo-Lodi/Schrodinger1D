package schrodinger.integrator;

import schrodinger.QuantumState;

public class RaptisAllison implements Integrator{
    /**
     * Propagates using the Raptis-Allison Exponentially Fitted Method.
     *
     * This method dynamically adjusts coefficients to be exact for the local frequency.
     * It behaves like Numerov when the energy is small/step is fine, but adapts
     * to high frequencies to maintain P-stability and accuracy.
     */
    public double propagate(double[] psi, int n, QuantumState state, Direction direction) {
        double h = state.getGrid().getStepSizeYCoordinate();
        double h2 = h * h;
        // 1. Identify indices
        int n_curr = n;
        int n_prev = n - direction.getValue();
        int n_next = n + direction.getValue();
        // 2. Calculate Q values (y'' = -Q y)
        double Q_next = state.QTildeValueAt(n_next);
        double Q_curr = state.QTildeValueAt(n_curr);
        double Q_prev = state.QTildeValueAt(n_prev);
        // 3. Calculate Z = h^2 * Q
        double Z_curr = h2 * Q_curr;
        // 4. Calculate Fitted Coefficients (beta and gamma)
        double beta, gamma;
        if (Math.abs(Z_curr) < 1e-4) {
            // Small Z limit: Standard Numerov Coefficients
            beta = 1.0 / 12.0;
            gamma = 10.0 / 12.0;
        } else if (Z_curr > 0) {
            // Oscillatory Region (Q > 0)
            double w = Math.sqrt(Z_curr);
            double halfW = w / 2.0;
            // Guard against tan(w/2) singularities (w near odd multiples of pi)
            if (Math.abs(Math.cos(halfW)) < 1e-10) {
                beta = 1.0 / 12.0;
                gamma = 10.0 / 12.0;
            } else {
                // Beta formula: 1/v^2 - 1/(2v tan(v/2))
                beta = (1.0 / Z_curr) - (1.0 / (2.0 * w * Math.tan(halfW)));
                // Gamma formula derived from consistency with exact solution y = exp(iw)
                // Solves: (2 - gamma*Z) = 2*cos(w)*(1 + beta*Z)
                gamma = (2.0 * Math.cos(w) * (1.0 + beta * Z_curr) - 2.0) / Z_curr;
            }
        } else {
            // Exponential Region (Q < 0)
            double w = Math.sqrt(-Z_curr);
            // Beta formula (switch to hyperbolic): 1/v^2 + 1/(2v tanh(v/2))
            // Note: Z_curr is negative, so 1/Z_curr is negative; the two terms partially cancel.
            beta = (1.0 / Z_curr) + (1.0 / (2.0 * w * Math.tanh(w / 2.0)));
            // Gamma formula for exponential (cosh instead of cos)
            gamma = (2.0 * Math.cosh(w) * (1.0 + beta * Z_curr) - 2.0) / Z_curr;
        }
        // 5. Propagate
        // Formula derived from y'' = -Qy substitution into the multistep form:
        // y_{n+1} (1 + beta*Z_next) = 2 y_n (1 + (gamma/2)*Z_n) - y_{n-1} (1 - beta*Z_{n-1})
        //
        // Note the minus sign in the previous term: (1 - beta*Z_prev).
        double numerator =
                2.0 * psi[n_curr] * (1.0 + 0.5 * gamma * Z_curr)
                        - psi[n_prev] * (1.0 + beta * (h2 * Q_prev));
        double denominator = 1.0 + beta * (h2 * Q_next);
        if (Math.abs(denominator) < 1e-14) {
            throw new ArithmeticException("Near-zero denominator in Raptis-Allison propagator at n=" + n);
        }
        return numerator / denominator;
    }


}
