package schrodinger.integrator.one_step;

import schrodinger.integrator.Integrator;

import java.util.function.DoubleUnaryOperator;

public class RK45DP implements Integrator {

    @Override
    public int minHistoryLength() {
        return 1;
    }

    @Override
    public int globalConvergenceOrder() {
        return 5;
    }


    //    It’s not “RKN” in the strict specialized sense, but it is tableau-based, one-step, and gives you 5th-order global accuracy with moderate code/coefficients.
    // standard 5th-order Runge–Kutta method (Dormand–Prince 5(4))
    // c nodes
    private final static double c2 = 1.0 / 5.0;
    private final static double c3 = 3.0 / 10.0;
    private final static double c4 = 4.0 / 5.0;
    private final static double c5 = 8.0 / 9.0;
    private final static double c6 = 1.0;
    private final static double c7 = 1.0;

    // a matrix (Dormand–Prince)
    private final static double a21 = 1.0 / 5.0;

    private final static double a31 = 3.0 / 40.0;
    private final static double a32 = 9.0 / 40.0;

    private final static double a41 = 44.0 / 45.0;
    private final static double a42 = -56.0 / 15.0;
    private final static double a43 = 32.0 / 9.0;

    private final static double a51 = 19372.0 / 6561.0;
    private final static double a52 = -25360.0 / 2187.0;
    private final static double a53 = 64448.0 / 6561.0;
    private final static double a54 = -212.0 / 729.0;

    private final static double a61 = 9017.0 / 3168.0;
    private final static double a62 = -355.0 / 33.0;
    private final static double a63 = 46732.0 / 5247.0;
    private final static double a64 = 49.0 / 176.0;
    private final static double a65 = -5103.0 / 18656.0;

    private final static double a71 = 35.0 / 384.0;
    private final static double a72 = 0.0;
    private final static double a73 = 500.0 / 1113.0;
    private final static double a74 = 125.0 / 192.0;
    private final static double a75 = -2187.0 / 6784.0;
    private final static double a76 = 11.0 / 84.0;

    // b weights for the 5th-order solution (same as a7 row)
    private final static double b1 = a71, b2 = a72, b3 = a73, b4 = a74, b5 = a75, b6 = a76, b7 = 0.0;

    @Override
    public double propagate(double[] psi, double[] currentPsiPrime, int n, double step,
                            DoubleUnaryOperator qTilde,
                            DoubleUnaryOperator qTildePrime,
                            DoubleUnaryOperator qTildeDoublePrime,
                            Direction direction) {

        final int d = direction.getValue();
        final double h = step * d;
        final int n0 = n;
        final int n1 = n + d;

        // "x" is your fractional index coordinate (since qTilde accepts fractional indices)
        final double x0 = n0;

        // state
        final double y0 = psi[n0];
        final double v0 = currentPsiPrime[0];

        // helper: acceleration v' = -q(x)*y
        // (q evaluated at fractional index x0 + ci*d, consistent with your convention)
        java.util.function.DoubleBinaryOperator acc = (x, y) -> -qTilde.applyAsDouble(x) * y;

        // k vectors: (ky_i, kv_i) = (y', v')
        final double ky1 = v0;
        final double kv1 = acc.applyAsDouble(x0, y0);

        final double y2 = y0 + h * (a21 * ky1);
        final double v2 = v0 + h * (a21 * kv1);
        final double ky2 = v2;
        final double kv2 = acc.applyAsDouble(x0 + c2 * d, y2);

        final double y3 = y0 + h * (a31 * ky1 + a32 * ky2);
        final double v3 = v0 + h * (a31 * kv1 + a32 * kv2);
        final double ky3 = v3;
        final double kv3 = acc.applyAsDouble(x0 + c3 * d, y3);

        final double y4 = y0 + h * (a41 * ky1 + a42 * ky2 + a43 * ky3);
        final double v4 = v0 + h * (a41 * kv1 + a42 * kv2 + a43 * kv3);
        final double ky4 = v4;
        final double kv4 = acc.applyAsDouble(x0 + c4 * d, y4);

        final double y5 = y0 + h * (a51 * ky1 + a52 * ky2 + a53 * ky3 + a54 * ky4);
        final double v5 = v0 + h * (a51 * kv1 + a52 * kv2 + a53 * kv3 + a54 * kv4);
        final double ky5 = v5;
        final double kv5 = acc.applyAsDouble(x0 + c5 * d, y5);

        final double y6 = y0 + h * (a61 * ky1 + a62 * ky2 + a63 * ky3 + a64 * ky4 + a65 * ky5);
        final double v6 = v0 + h * (a61 * kv1 + a62 * kv2 + a63 * kv3 + a64 * kv4 + a65 * kv5);
        final double ky6 = v6;
        final double kv6 = acc.applyAsDouble(x0 + c6 * d, y6);

        final double y7 = y0 + h * (a71 * ky1 + a73 * ky3 + a74 * ky4 + a75 * ky5 + a76 * ky6);
        final double v7 = v0 + h * (a71 * kv1 + a73 * kv3 + a74 * kv4 + a75 * kv5 + a76 * kv6);
        final double ky7 = v7;
        final double kv7 = acc.applyAsDouble(x0 + c7 * d, y7);

        // 5th-order update
        final double yNext = y0 + h * (b1 * ky1 + b2 * ky2 + b3 * ky3 + b4 * ky4 + b5 * ky5 + b6 * ky6 + b7 * ky7);
        final double vNext = v0 + h * (b1 * kv1 + b2 * kv2 + b3 * kv3 + b4 * kv4 + b5 * kv5 + b6 * kv6 + b7 * kv7);

        currentPsiPrime[0] = vNext;
        return yNext;
    }

    @Override
    public double[] getFractionalOffsets() {
        // it also calls 1/9 (=1-c5) and 7/10 (1-c3) when going backwards
        return new double[]{0.0, 1. - c5, c2, c3, 1. - c3, c4, c5};
    }

}

