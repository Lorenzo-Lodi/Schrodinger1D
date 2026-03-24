package schrodinger.integrator.rungekutta;

import schrodinger.integrator.Integrator;

import java.util.function.DoubleUnaryOperator;

public final class CFMagnus6 implements Integrator {

    // 6th order requires 3 Gauss-Legendre nodes and 4 exponentials
    private static final int NODES = 3;
    private static final int EXPONENTIALS = 4;

    private static final double[] C = new double[NODES];
    private static final double[] GW = new double[NODES];

    // F arrays: 4 exponentials, 3 Legendre coefficients (A_1, A_2, A_3)
    private static final double[][] F = new double[EXPONENTIALS][3];

    private static final double[] A1W = new double[EXPONENTIALS];
    private static final double[][] V = new double[EXPONENTIALS][NODES];

    static {
        // 1. Gauss-Legendre Nodes mapped to [0, 1]
        double sqrt35 = Math.sqrt(3.0 / 5.0);
        C[0] = 0.5 - 0.5 * sqrt35;
        C[1] = 0.5;
        C[2] = 0.5 + 0.5 * sqrt35;

        // 2. Quadrature weights for interval size 1
        GW[0] = 5.0 / 18.0;
        GW[1] = 8.0 / 18.0;  // or 4.0 / 9.0
        GW[2] = 5.0 / 18.0;

        // 3. User-provided coefficients for CF6:4 (First half)
        // From Table 3 of 1102.5071v2.pdf
        F[0][0] = 1.0798524263824308825;
        F[0][1] = F[0][0] - (2.0 / 3.0) * F[0][0] * F[0][0];
        F[0][2] = 1.0 / (10.0 - 10.0 * F[0][0]);

        F[1][0] = 0.5 - F[0][0];
        F[1][1] = (1.0 - 4.0 * F[0][0] + 2.0 * F[0][0] * F[0][0]) / 3.0;
        F[1][2] = -F[0][2];

        // 4. Apply Time-Reversal Symmetry for the second half
        // F[2] is the mirror of F[1], F[3] is the mirror of F[0]
        for (int i = 0; i < 2; i++) {
            F[3 - i][0] = F[i][0]; // Even parity
            F[3 - i][1] = -F[i][1]; // Odd parity flips
            F[3 - i][2] = F[i][2]; // Even parity
        }

        // 5. Precompute effective node weights using exact polynomials
        for (int i = 0; i < EXPONENTIALS; i++) {
            A1W[i] = F[i][0];

            for (int k = 0; k < NODES; k++) {
                double c = C[k];

                double p0 = 1.0;
                double p1 = 2.0 * c - 1.0;
                double p2 = 6.0 * c * c - 6.0 * c + 1.0;

                V[i][k] = GW[k] * (
                        F[i][0] * p0
                                + 3.0 * F[i][1] * p1
                                + 5.0 * F[i][2] * p2
                );
            }
        }

    }

    @Override
    public int minHistoryLength() {
        return 1;
    }

    @Override
    public int globalConvergenceOrder() {
        return 6;
    }

    @Override
    public double propagate(double[] psi, double[] currentPsiPrime, int n, double step,
                            DoubleUnaryOperator qTilde,
                            DoubleUnaryOperator qTildePrime,
                            DoubleUnaryOperator qTildeDoublePrime,
                            Direction direction) {

        double hd = step * direction.getValue();
        double y = psi[n];
        double yp = currentPsiPrime[0];

        // 1. Evaluate Q(x) at the Gauss-Legendre nodes
        double[] Q = new double[NODES];
        for (int k = 0; k < NODES; k++) {
            Q[k] = qTilde.applyAsDouble(n + C[k] * direction.getValue());
        }

        // 2. Multiply the state vector by the sequence of Exponentials (Must be in reverse order!)
        for (int i = EXPONENTIALS - 1; i >= 0; i--) {
            double qComb = 0.0;
            for (int k = 0; k < NODES; k++) {
                qComb += V[i][k] * Q[k];
            }

            double a = hd * A1W[i];
            double b = hd * qComb;

            double z = a * b;
            double yNext, ypNext;

            if (z > 0.0) {
                // Oscillatory regime (Q > 0)
                double w = Math.sqrt(z);
                double cos = Math.cos(w);
                double sinc = Math.sin(w) / w;

                yNext = cos * y + a * sinc * yp;
                ypNext = -b * sinc * y + cos * yp;
            } else if (z < 0.0) {
                // Exponential regime (Q < 0)
                double nu = Math.sqrt(-z);
                double cosh = Math.cosh(nu);
                double sinhc = Math.sinh(nu) / nu;

                yNext = cosh * y + a * sinhc * yp;
                ypNext = -b * sinhc * y + cosh * yp;
            } else {
                // Edge case: Q = 0 (Free particle)
                yNext = y + a * yp;
                ypNext = yp - b * y;
            }

            y = yNext;
            yp = ypNext;
        }

        currentPsiPrime[0] = yp;
        return y;
    }
}
