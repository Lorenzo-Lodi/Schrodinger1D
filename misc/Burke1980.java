package schrodinger.integrator;

import java.util.function.IntToDoubleFunction;

public class Burke1980 implements Integrator {
// This is formally correct but the underlying method is unstable (Dahlquist's Root Condition) and unusable

    @Override
    public int minHistoryLength() {
        return 4;
    }

    @Override
    public int globalConvergenceOrder() {
        return 4;
    }

    @Override
    public double propagate(double[] psi, int n, double step, IntToDoubleFunction qTilde, Direction direction) {
        double h2 = step * step;
        int d = direction.getValue();

        // Target point to calculate
        int nP1 = n + d;

        // Past known points
        int n0 = n;
        int n1 = n - d;
        int n2 = n - 2 * d;
        int n3 = n - 3 * d;

        final double Q_nP1 = qTilde.applyAsDouble(nP1);
        final double Q_n0  = qTilde.applyAsDouble(n0);
        final double Q_n1  = qTilde.applyAsDouble(n1);
        final double Q_n2  = qTilde.applyAsDouble(n2);
        final double Q_n3  = qTilde.applyAsDouble(n3);

        double num = -465. * psi[n3] - 1920. * psi[n2] + 4770. * psi[n1] - 1920. * psi[n0]
                - h2 * (688. * Q_n0 * psi[n0] + 2358. * Q_n1 * psi[n1] + 688. * Q_n2 * psi[n2] + 23. * Q_n3 * psi[n3]);

        double denom = 465. + 23. * h2 * Q_nP1;

        return num / denom;
    }


}
