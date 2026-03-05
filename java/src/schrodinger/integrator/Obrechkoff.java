package schrodinger.integrator;

import java.util.function.DoubleUnaryOperator;

public class Obrechkoff implements Integrator {

    @Override
    public int minHistoryLength() {
        return 2;
    }

    @Override
    public int globalConvergenceOrder() {
        return 6;
    }


    // https://www.perplexity.ai/search/can-an-adams-moulton-method-be-4iA7wTxyTHOyelu30RUUGA#2254f9eb-2c6a-49ad-86e3-c5fbae55a1b0
    @Override
    public double propagate(double[] psi, int n, double step, DoubleUnaryOperator qTilde, Direction direction) {
        int d = direction.getValue();
        // TODO !!!!
        double Q_prime_n = 0;
        double Q_prime_prime_n = 0;
        // TODO !!!!

        double Q_n = qTilde.applyAsDouble(n);
        double Q_np1 = qTilde.applyAsDouble(n + d);
        double Q_nm1 = qTilde.applyAsDouble(n - d);

        double h2 = step * step;
        double h3 = h2 * step;
        double C_n = 3. * h3 * Q_prime_n / (10. * (6. - h2 * Q_n));
        double A_n = 1. + (h2 / 30.) * Q_np1 + C_n;
        double B_n = 2. - (14. / 15.) * h2 * Q_n - C_n * h3 * Q_prime_n / 3. + (h2 * h2 / 20.) * (Q_n * Q_n - Q_prime_prime_n);
        double E_n = 1. + (h2 / 30.) * Q_nm1 - C_n;
        return (B_n * psi[n] - E_n * psi[n - d]) / A_n;
    }

}
