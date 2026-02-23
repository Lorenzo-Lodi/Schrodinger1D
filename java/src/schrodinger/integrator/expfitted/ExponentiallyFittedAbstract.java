package schrodinger.integrator.expfitted;

import schrodinger.integrator.Integrator;

import java.util.function.IntToDoubleFunction;

/**
 * Abstract base class for exponentially fitted integration methods.
 */
public abstract class ExponentiallyFittedAbstract implements Integrator {

    @Override
    public double propagate(double[] psi, int n, double step, IntToDoubleFunction qTilde, Direction direction) {
        double h = step;
        double h2 = h * h;

        // 1. Identify neighbors
        int prev = n - direction.getValue();
        int curr = n;
        int next = n + direction.getValue();

        // 2. Calculate Z = h^2 * Q for each point
        // Positive Q -> Oscillatory. Negative Q -> Exponential.
        double Z_prev = h2 * qTilde.applyAsDouble(prev);
        double Z_curr = h2 * qTilde.applyAsDouble(curr);
        double Z_next = h2 * qTilde.applyAsDouble(next);

        // 3. Calculate EF coefficients dynamically
        double beta_prev = getBeta(Z_prev);
        double beta_next = getBeta(Z_next);
        double beta_curr = getBeta(Z_curr);
        double gamma_curr = getGamma(Z_curr, beta_curr);

        // 4. Propagate: (1 + β Z_next) y_next = (2 - γ Z_curr) y_curr - (1 + β Z_prev) y_prev
        double numerator = (2.0 - gamma_curr * Z_curr) * psi[curr] - (1.0 + beta_prev * Z_prev) * psi[prev];
        double denominator = 1.0 + beta_next * Z_next;

        return numerator / denominator;
    }

    @Override
    public int minHistoryLength() { return 2; }

    @Override
    public int globalConvergenceOrder() { return 4; }

    /**
     * Returns the beta coefficient for exponentially fitted method given Z = h^2 * Q.
     */
    public abstract double getBeta(double Z);

    /**
     * Returns the gamma coefficient for exponentially fitted method given Z = h^2 * Q and beta.
     */
    public abstract double getGamma(double Z, double beta);
}
