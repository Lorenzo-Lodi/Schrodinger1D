package schrodinger.integrator.rungekutta;

import schrodinger.integrator.Integrator;

import java.util.function.DoubleUnaryOperator;

public final class CFMagnus8 implements Integrator {

    private static final int NODES = 4;
    private static final int EXPONENTIALS = 11;

    private static final double[] C = new double[NODES]; // Gauss-Legendre nodes
    private static final double[][] W = new double[EXPONENTIALS][NODES]; // Weights

    static {
        // 1. Compute the 4 Gauss-Legendre Nodes on [0, 1]
        double sqrt12 = Math.sqrt(1.2);
        double term1 = 0.5 * Math.sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt12);
        double term2 = 0.5 * Math.sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt12);
        C[0] = 0.5 - term1;
        C[1] = 0.5 - term2;
        C[2] = 0.5 + term2;
        C[3] = 0.5 + term1;

        // i = 1 (First exponential)
        W[0][0] = 0.169715531043933180094151;
        W[0][1] = 0.152866146944615909929839;
        W[0][2] = 0.119167378745981369601216;
        W[0][3] = 0.068619226448029559107538;

        // i = 2
        W[1][0] = 0.379420807516005431504230;
        W[1][1] = 0.148839980923180990943008;
        W[1][2] = -0.115880829186628075021088;
        W[1][3] = -0.188555246668412628269760;

        // i = 3
        W[2][0] = 0.469459306644050573017994;
        W[2][1] = -0.379844237839363505173921;
        W[2][2] = 0.022898814729462898505141;
        W[2][3] = 0.571855043580130805495594;

        // i = 4
        W[3][0] = -0.448225927391070886302766;
        W[3][1] = 0.362889857410989942809900;
        W[3][2] = -0.022565582830528472333301;
        W[3][3] = -0.544507517141613383517695;

        // i = 5
        W[4][0] = -0.293924473106317605373923;
        W[4][1] = -0.026255628265819381983204;
        W[4][2] = 0.096761509131620390100068;
        W[4][3] = 0.000018330145571671744069;

        // i = 6 (Center exponential)
        W[5][0] = 0.447109510586798614120629;
        W[5][2] = -0.200762581179816221704073;
        W[5][1] = W[5][2];
        W[5][3] = W[5][0];

        // 3. Automatically mirror the coefficients for the second half (i = 7 to 11)
        // Time-reversal symmetry dictates: f_{12-i, 5-j} = f_{i, j}
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 4; j++) {
                W[10 - i][3 - j] = W[i][j];
            }
        }

    }

    @Override
    public int minHistoryLength() {
        return 1;
    }

    @Override
    public int globalConvergenceOrder() {
        return 8;

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
