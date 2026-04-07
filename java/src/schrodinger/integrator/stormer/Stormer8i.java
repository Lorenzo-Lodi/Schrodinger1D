package schrodinger.integrator.stormer;

import schrodinger.integrator.Integrator;

import java.util.function.DoubleUnaryOperator;

public class Stormer8i implements Integrator {

    @Override
    public int minHistoryLength() {
        return 7;
    }

    @Override
    public int globalConvergenceOrder() {
        return 8;
    }

    @Override
    public boolean needsPotentialCapping() {
        return true;
    }

    // 8th-order implicit Cowell coefficients
    private final static double c_next = 275.0 / 4032.0;
    private final static double c_curr = 13831.0 / 15120.0;
    private final static double c_prev = -2099.0 / 20160.0;
    private final static double c_prev2 = 811.0 / 3360.0;
    private final static double c_prev3 = -11477.0 / 60480.0;
    private final static double c_prev4 = 29.0 / 315.0;
    private final static double c_prev5 = -517.0 / 20160.0;
    private final static double c_prev6 = 19.0 / 6048.0;

    @Override
    public double propagate(double[] psi, double[] currentPsiPrime, int n, double step,
                            DoubleUnaryOperator qTilde,
                            DoubleUnaryOperator qTildePrime,
                            DoubleUnaryOperator qTildeDoublePrime, Direction direction) {

        double h2 = step * step;
        int d = direction.getValue();

        int n_next = n + d;
        int n0 = n;
        int n1 = n - d;
        int n2 = n - 2 * d;
        int n3 = n - 3 * d;
        int n4 = n - 4 * d;
        int n5 = n - 5 * d;
        int n6 = n - 6 * d;

        // Evaluate qTilde at the future point (implicit solved algebraically)
        double Q_next = qTilde.applyAsDouble(n_next);
        double Q_curr = qTilde.applyAsDouble(n0);
        double Q_prev = qTilde.applyAsDouble(n1);
        double Q_prev2 = qTilde.applyAsDouble(n2);
        double Q_prev3 = qTilde.applyAsDouble(n3);
        double Q_prev4 = qTilde.applyAsDouble(n4);
        double Q_prev5 = qTilde.applyAsDouble(n5);
        double Q_prev6 = qTilde.applyAsDouble(n6);

        // Explicit sum of the history points
        double rhs = 2.0 * psi[n0] - psi[n1]
                - h2 * (c_curr * Q_curr * psi[n0]
                + c_prev * Q_prev * psi[n1]
                + c_prev2 * Q_prev2 * psi[n2]
                + c_prev3 * Q_prev3 * psi[n3]
                + c_prev4 * Q_prev4 * psi[n4]
                + c_prev5 * Q_prev5 * psi[n5]
                + c_prev6 * Q_prev6 * psi[n6]);

        // Divide by the algebraic future term to complete the implicit step
        return rhs / (1.0 + h2 * c_next * Q_next);
    }
}
