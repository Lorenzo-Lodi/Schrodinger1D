package schrodinger.integrator.misc;

import schrodinger.integrator.Integrator;

import java.util.function.DoubleUnaryOperator;

public class Obrechkoff6 implements Integrator {

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
    public double propagate(double[] psi, double[] currentPsiPrime, int n, double step,
                            DoubleUnaryOperator qTilde,
                            DoubleUnaryOperator qTildePrime,
                            DoubleUnaryOperator qTildeDoublePrime, Direction direction) {
        int d = direction.getValue();

        double Q_prime_n = qTildePrime.applyAsDouble(n);
        double Q_prime_prime_n = qTildeDoublePrime.applyAsDouble(n);

        double Q_n = qTilde.applyAsDouble(n);
        double Q_np1 = qTilde.applyAsDouble(n + d);
        double Q_nm1 = qTilde.applyAsDouble(n - d);

        // 1. Define H as signed step
        double H = d * step;
        double H2 = H * H;
        double H3 = H2 * H;

        double C_n = 3. * H3 * Q_prime_n / (10. * (6. - H2 * Q_n));
        double A_n = 1. + (H2 / 30.) * Q_np1 + C_n;
        double E_n = 1. + (H2 / 30.) * Q_nm1 - C_n;
        double B_n = 2. - (14. / 15.) * H2 * Q_n - C_n * H3 * Q_prime_n / 3.
                + (H2 * H2 / 20.) * (Q_n * Q_n - Q_prime_prime_n);

        return (B_n * psi[n] - E_n * psi[n - d]) / A_n;

    }

}
