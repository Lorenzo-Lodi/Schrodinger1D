package schrodinger.integrator;

import schrodinger.QuantumState;

public abstract class ExponentiallyFittedAbstract implements Integrator {

    public double propagate(double[] psi, int n, QuantumState state, Direction direction) {
        double h = state.getGrid().getStepSizeYCoordinate();
        double h2 = h * h;

        // 1. Identify neighbors
        int prev = n - direction.getValue();
        int curr = n;
        int next = n + direction.getValue();

        // 2. Calculate Z = h^2 * Q for each point
        // Positive Q -> Oscillatory. Negative Q -> Exponential.
        double Z_prev = h2 * state.QTildeValueAt(prev);
        double Z_curr = h2 * state.QTildeValueAt(curr);
        double Z_next = h2 * state.QTildeValueAt(next);

        // 3. Calculate EF coefficients dynamically
        double beta_prev = getBeta(Z_prev);
        double beta_next = getBeta(Z_next);
        double beta_curr = getBeta(Z_curr);
        double gamma_curr = getGamma(Z_curr, beta_curr);

        // 4. Propagate: (1 + β Z_next) y_next = (2 - γ Z_curr) y_curr - (1 + β Z_prev) y_prev
        double numerator = (2.0 - Z_curr * gamma_curr) * psi[curr]
                - (1.0 + Z_prev * beta_prev) * psi[prev];

        double denominator = 1.0 + Z_next * beta_next;

        return numerator / denominator;
    }

    public abstract double getBeta(double Z);

    public abstract double getGamma(double Z, double beta);
}
