package schrodinger.integrator.rungekutta;

import schrodinger.integrator.Integrator;

import java.util.function.DoubleUnaryOperator;

public final class CFMagnus4 implements Integrator {

    private static final int NODES = 2;
    private static final int EXPONENTIALS = 2;

    private static final double[] C = new double[NODES]; // Gauss-Legendre nodes
    private static final double[][] W = new double[EXPONENTIALS][NODES]; // Weights

    static {
        double sqrt3 = Math.sqrt(3.0);

        // 1. Gauss-Legendre Nodes on [0, 1]
        C[0] = 0.5 - sqrt3 / 6.0;
        C[1] = 0.5 + sqrt3 / 6.0;

        // 2. Weights for Exponential 1
        W[0][0] = 0.25 + sqrt3 / 6.0;
        W[0][1] = 0.25 - sqrt3 / 6.0;

        // 3. Weights for Exponential 2
        W[1][0] = 0.25 - sqrt3 / 6.0;
        W[1][1] = 0.25 + sqrt3 / 6.0;
    }

    @Override
    public int minHistoryLength() {
        return 1;
    }

    @Override
    public int globalConvergenceOrder() {
        return 4; // Change to 8 when using CF8 coefficients
    }

    @Override
    public double propagate(double[] psi, double[] currentPsiPrime, int n, double step,
                            DoubleUnaryOperator qTilde,
                            DoubleUnaryOperator qTildePrime,
                            DoubleUnaryOperator qTildeDoublePrime,
                            Direction direction) {

        double h = step;
        int d = direction.getValue();

        double y = psi[n];
        double yp = currentPsiPrime[0];

        // 1. Pre-calculate Q(x) at the Gauss-Legendre quadrature nodes
        double[] Q = new double[NODES];
        for (int k = 0; k < NODES; k++) {
            Q[k] = qTilde.applyAsDouble(n + C[k] * d);
        }

        // 2. Multiply the state vector by the sequence of Exponentials
        for (int i = 0; i < EXPONENTIALS; i++) {

            double sumW = 0.0;
            double sumWQ = 0.0;

            for (int k = 0; k < NODES; k++) {
                sumW += W[i][k];
                sumWQ += W[i][k] * Q[k];
            }

            // Elements of the combined matrix C_i
            double W_i = h * d * sumW;
            double K_i = h * d * sumWQ;

            // 3. Exact Analytical 2x2 Matrix Exponential
            double omegaSq = W_i * K_i;
            double yNext, ypNext;

            if (omegaSq > 0) {
                // Oscillatory regime (Q > 0)
                double omega = Math.sqrt(omegaSq);
                double cos = Math.cos(omega);
                double sinc = Math.sin(omega) / omega; // sin(w)/w

                yNext = cos * y + W_i * sinc * yp;
                ypNext = -K_i * sinc * y + cos * yp;
            } else if (omegaSq < 0) {
                // Exponential regime (Q < 0)
                double nu = Math.sqrt(-omegaSq);
                double cosh = Math.cosh(nu);
                double sinhc = Math.sinh(nu) / nu; // sinh(v)/v

                yNext = cosh * y + W_i * sinhc * yp;
                ypNext = -K_i * sinhc * y + cosh * yp;
            } else {
                // Edge case: Q = 0 (Free particle)
                yNext = y + W_i * yp;
                ypNext = yp - K_i * y;
            }

            y = yNext;
            yp = ypNext;
        }

        currentPsiPrime[0] = yp;
        return y;
    }
}
