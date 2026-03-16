package schrodinger.integrator.rungekutta;

import schrodinger.integrator.Integrator;

import java.util.function.DoubleUnaryOperator;

public final class RKN4 implements Integrator {

    private static final double C2 = 0.5;
    private static final double C3 = 1.0;

    private static final double A21 = 1.0 / 8.0;
    private static final double A31 = 0.0;
    private static final double A32 = 1.0 / 2.0;

    // weights for y update (position)
    private static final double B1 = 1.0 / 6.0;
    private static final double B2 = 1.0 / 3.0;
    private static final double B3 = 0.0;

    // weights for y' update (velocity)
    private static final double BP1 = 1.0 / 6.0;
    private static final double BP2 = 2.0 / 3.0;
    private static final double BP3 = 1.0 / 6.0;

    @Override
    public int minHistoryLength() {
        return 1;
    }

    @Override
    public int globalConvergenceOrder() {
        return 4;
    }

    @Override
    public double propagate(double[] psi, double[] currentPsiPrime, int n, double step,
                            DoubleUnaryOperator qTilde,
                            DoubleUnaryOperator qTildePrime,
                            DoubleUnaryOperator qTildeDoublePrime,
                            Direction direction) {

        double h = step;
        double h2 = step * step;
        int d = direction.getValue();

        int n0 = n;
        int n1 = n + d; // target index (n+1)

        double q1 = qTilde.applyAsDouble(n0);
        double f1 = -q1 * psi[n0];

        double Y2 = psi[n0] + C2 * h * d * currentPsiPrime[0] + h2 * A21 * f1;
        double q2 = qTilde.applyAsDouble(n0 + C2 * d);
        double f2 = -q2 * Y2;

        double Y3 = psi[n0] + C3 * h * d * currentPsiPrime[0] + h2 * (A31 * f1 + A32 * f2);
        double q3 = qTilde.applyAsDouble(n0 + C3 * d);
        double f3 = -q3 * Y3;

        // advance y and y'
        double result = psi[n0] + h * d * currentPsiPrime[0] + h2 * (B1 * f1 + B2 * f2 + B3 * f3);
        currentPsiPrime[0] = currentPsiPrime[0] + h * d * (BP1 * f1 + BP2 * f2 + BP3 * f3);
        return result;
    }
}
