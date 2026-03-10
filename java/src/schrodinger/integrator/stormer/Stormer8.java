package schrodinger.integrator.stormer;

import schrodinger.integrator.Integrator;

import java.util.function.DoubleUnaryOperator;

public class Stormer8 implements Integrator {

    @Override
    public int minHistoryLength() {
        return 8;
    }

    @Override
    public int globalConvergenceOrder() {
        return 8;
    }

    private final static double b_curr = 22081.0 / 15120.0;
    private final static double b_prev = 4511.0 / 2240.0;
    private final static double b_prev2 = 40933.0 / 10080.0;
    private final static double b_prev3 = 300227.0 / 60480.0;
    private final static double b_prev4 = 9857.0 / 2520.0;
    private final static double b_prev5 = 39017.0 / 20160.0;
    private final static double b_prev6 = 3319.0 / 6048.0;
    private final static double b_prev7 = 275.0 / 4032;

    @Override
    public double propagate(double[] psi, double[] currentPsiPrime, int n, double step,
                            DoubleUnaryOperator qTilde,
                            DoubleUnaryOperator qTildePrime,
                            DoubleUnaryOperator qTildeDoublePrime, Direction direction) {
        double h2 = step * step;
        int d = direction.getValue();

        int n0 = n;
        int n1 = n - d;
        int n2 = n - 2 * d;
        int n3 = n - 3 * d;
        int n4 = n - 4 * d;
        int n5 = n - 5 * d;
        int n6 = n - 6 * d;
        int n7 = n - 7 * d;

        double Q_curr = qTilde.applyAsDouble(n0);
        double Q_prev = qTilde.applyAsDouble(n1);
        double Q_prev2 = qTilde.applyAsDouble(n2);
        double Q_prev3 = qTilde.applyAsDouble(n3);
        double Q_prev4 = qTilde.applyAsDouble(n4);
        double Q_prev5 = qTilde.applyAsDouble(n5);
        double Q_prev6 = qTilde.applyAsDouble(n6);
        double Q_prev7 = qTilde.applyAsDouble(n7);

        return 2.0 * psi[n0] - psi[n1]
                - h2 * (b_curr * Q_curr * psi[n0]
                - b_prev * Q_prev * psi[n1]
                + b_prev2 * Q_prev2 * psi[n2]
                - b_prev3 * Q_prev3 * psi[n3]
                + b_prev4 * Q_prev4 * psi[n4]
                - b_prev5 * Q_prev5 * psi[n5]
                + b_prev6 * Q_prev6 * psi[n6]
                - b_prev7 * Q_prev7 * psi[n7]);
    }


}
