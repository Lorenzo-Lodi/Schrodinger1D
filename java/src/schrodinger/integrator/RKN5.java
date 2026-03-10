package schrodinger.integrator;

import java.util.function.DoubleUnaryOperator;

public class RKN5 implements Integrator {

    // Stage nodes (c_i)
    private static final double C2 = 0.5; // Actually 1/2 but written nicely
    private static final double C3 = 0.5 + Math.sqrt(5.0) / 10.0; // ~0.7236
    private static final double C4 = 0.5 - Math.sqrt(5.0) / 10.0; // ~0.2763

    // Weights for updating Position (b_i)
    private static final double B1 = 1.0 / 12.0;
    private static final double B2 = 5.0 / 12.0;
    private static final double B3 = 5.0 / 12.0;
    private static final double B4 = 1.0 / 12.0;

    // Weights for updating Velocity (b'_i)
    private static final double BP1 = 1.0 / 12.0;
    private static final double BP2 = 5.0 / 12.0;
    private static final double BP3 = 5.0 / 12.0;
    private static final double BP4 = 1.0 / 12.0;

    // Internal matrix coefficients (a_ij)
    private static final double A21 = 1.0 / 8.0;
    private static final double A31 = (1.0 + Math.sqrt(5.0)) / 10.0;
    private static final double A32 = (1.0 - Math.sqrt(5.0)) / 10.0;
    private static final double A41 = (1.0 - Math.sqrt(5.0)) / 10.0;
    private static final double A42 = (1.0 + Math.sqrt(5.0)) / 10.0;
    private static final double A43 = 0.0;

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
    @Override
    public double propagate(double[] psi, double[] psiPrime, int n, double step,
                            DoubleUnaryOperator qTilde,
                            DoubleUnaryOperator qTildePrime,
                            DoubleUnaryOperator qTildeDoublePrime,
                            Direction direction) {

        final double h = step;
        final int d = direction.getValue();
        final int n0 = n;
        final int n1 = n + d;

        // "x" is your fractional index coordinate (since qTilde accepts fractional indices)
        final double x0 = n0;

        // c nodes
        final double c2 = 1.0 / 5.0;
        final double c3 = 3.0 / 10.0;
        final double c4 = 4.0 / 5.0;
        final double c5 = 8.0 / 9.0;
        final double c6 = 1.0;
        final double c7 = 1.0;

        // a matrix (Dormand–Prince)
        final double a21 = 1.0 / 5.0;

        final double a31 = 3.0 / 40.0;
        final double a32 = 9.0 / 40.0;

        final double a41 = 44.0 / 45.0;
        final double a42 = -56.0 / 15.0;
        final double a43 = 32.0 / 9.0;

        final double a51 = 19372.0 / 6561.0;
        final double a52 = -25360.0 / 2187.0;
        final double a53 = 64448.0 / 6561.0;
        final double a54 = -212.0 / 729.0;

        final double a61 = 9017.0 / 3168.0;
        final double a62 = -355.0 / 33.0;
        final double a63 = 46732.0 / 5247.0;
        final double a64 = 49.0 / 176.0;
        final double a65 = -5103.0 / 18656.0;

        final double a71 = 35.0 / 384.0;
        final double a72 = 0.0;
        final double a73 = 500.0 / 1113.0;
        final double a74 = 125.0 / 192.0;
        final double a75 = -2187.0 / 6784.0;
        final double a76 = 11.0 / 84.0;

        // b weights for the 5th-order solution (same as a7 row)
        final double b1 = a71, b2 = a72, b3 = a73, b4 = a74, b5 = a75, b6 = a76, b7 = 0.0;

        // state
        final double y0 = psi[n0];
        final double v0 = psiPrime[n0];

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

        psiPrime[n1] = vNext;
        return yNext;
    }
}

