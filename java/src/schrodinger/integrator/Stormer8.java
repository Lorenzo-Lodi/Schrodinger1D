package schrodinger.integrator;

import java.util.function.IntToDoubleFunction;

public class Stormer8 implements Integrator {

    @Override
    public int minHistoryLength() {
        return 8;
    }

    @Override
    public int globalConvergenceOrder() {
        return 8;
    }

    @Override
    public double propagate(double[] psi, int n, double step, IntToDoubleFunction qTilde, Direction direction) {
        double h2 = step * step;
        int d = direction.getValue();

        // Indices of the 8 past points
        int n0 = n, n1 = n - d, n2 = n - 2 * d, n3 = n - 3 * d;
        int n4 = n - 4 * d, n5 = n - 5 * d, n6 = n - 6 * d, n7 = n - 7 * d;

        // Q evaluated at integer grid points (no fractional positions needed)
        double F0 = qTilde.applyAsDouble(n0) * psi[n0];
        double F1 = qTilde.applyAsDouble(n1) * psi[n1];
        double F2 = qTilde.applyAsDouble(n2) * psi[n2];
        double F3 = qTilde.applyAsDouble(n3) * psi[n3];
        double F4 = qTilde.applyAsDouble(n4) * psi[n4];
        double F5 = qTilde.applyAsDouble(n5) * psi[n5];
        double F6 = qTilde.applyAsDouble(n6) * psi[n6];
        double F7 = qTilde.applyAsDouble(n7) * psi[n7];

        return 2.0 * psi[n0] - psi[n1]
                - h2 * ((22081.0 / 15120) * F0
                - (4511.0 / 2240) * F1
                + (40933.0 / 10080) * F2
                - (300227.0 / 60480) * F3
                + (9857.0 / 2520) * F4
                - (39017.0 / 20160) * F5
                + (3319.0 / 6048) * F6
                - (275.0 / 4032) * F7);
    }


}
