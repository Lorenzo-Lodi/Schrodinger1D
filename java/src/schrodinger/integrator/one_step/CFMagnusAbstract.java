package schrodinger.integrator.one_step;

import schrodinger.integrator.Integrator;

import java.util.function.DoubleUnaryOperator;

public abstract class CFMagnusAbstract {

    private final int NODES;
    private final int EXPONENTIALS;
    private final double[] C;
    private final double[][] W;
    private final double[] A1W;
    private final double[] Q; // pre-allocated to avoid per-step allocation


    public CFMagnusAbstract(int nodes, int exponentials, double[] c, double[][] w, double[] a1W) {
        NODES = nodes;
        EXPONENTIALS = exponentials;
        C = c;
        W = w;
        A1W = a1W;
        Q = new double[NODES];
    }

    public double propagate(double[] psi, double[] currentPsiPrime, int n, double step,
                            DoubleUnaryOperator qTilde,
                            DoubleUnaryOperator qTildePrime,
                            DoubleUnaryOperator qTildeDoublePrime,
                            Integrator.Direction direction) {

        double hd = step * direction.getValue();
        double y = psi[n];
        double yp = currentPsiPrime[0];

        for (int k = 0; k < NODES; k++) {
            Q[k] = qTilde.applyAsDouble(n + C[k] * direction.getValue());
        }

        // 2. Multiply the state vector by the sequence of Exponentials (Must be in reverse order!)
        for (int i = EXPONENTIALS - 1; i >= 0; i--) {

            double sumWQ = 0.0;
            for (int k = 0; k < NODES; k++) {
                sumWQ += W[i][k] * Q[k];
            }

            double W_i = hd * A1W[i];
            double K_i = hd * sumWQ;

            double z = W_i * K_i;
            double yNext, ypNext;

            if (z > 0.0) {
                // Oscillatory regime (Q > 0)
                double w = Math.sqrt(z);
                double cos = Math.cos(w);
                double sinc = Math.sin(w) / w;

                yNext = cos * y + W_i * sinc * yp;
                ypNext = -K_i * sinc * y + cos * yp;
            } else if (z < 0.0) {
                // Exponential regime (Q < 0)
                double nu = Math.sqrt(-z);
                double cosh = Math.cosh(nu);
                double sinhc = Math.sinh(nu) / nu;

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


    /**
     * Precomputes effective node weights using shifted Legendre polynomials on [0,1], and
     * also copies the first coefficient of the i‑th exponential into the A1W vector.
     * <p>
     * For each exponential index {@code i} and each node {@code k}, this method calculates
     * the weight matrix entry {@code W[i][k]} as the product of the Gauss weight {@code GW[k]}
     * and a weighted sum of shifted Legendre polynomials evaluated at {@code x = C[k]}.
     * The coefficients of the sum are taken from {@code F[i][n]} and multiplied by the
     * odd factors {@code (2n+1)}. The number of polynomial terms used is determined by
     * the length of {@code F[i]} (i.e., {@code maxOrder = F[i].length - 1}),
     * allowing the same method to handle both the 3‑term and 4‑term variants.
     *
     * @param F            2‑D array where {@code F[i][n]} provides the coefficient for the
     *                     {@code n}-th shifted Legendre polynomial in the {@code i}-th exponential.
     * @param C            Array of node coordinates (values in [0,1]).
     * @param GW           Array of Gauss weights corresponding to each node.
     * @param W            Output matrix to be filled with the computed weights.
     * @param A1W          Output vector receiving {@code F[i][0]} for each {@code i}.
     * @param EXPONENTIALS Number of exponential terms (outer loop size).
     * @param NODES        Number of nodes (inner loop size).
     */
    public static void computeW(double[][] F, double[] C, double[] GW,
                                double[][] W, double[] A1W,
                                int EXPONENTIALS, int NODES) {

        for (int i = 0; i < EXPONENTIALS; i++) {
            A1W[i] = F[i][0];                     // common first step

            // Determine how many terms we have for this i (based on F[i] length)
            int maxOrder = F[i].length - 1;       // because F[i][0] corresponds to n=0

            for (int k = 0; k < NODES; k++) {
                double x = C[k];                  // renamed from c/x for clarity

                // Compute shifted Legendre polynomials up to maxOrder
                double[] P = new double[maxOrder + 1];
                if (maxOrder >= 0) P[0] = 1.0;                     // P0
                if (maxOrder >= 1) P[1] = 2.0 * x - 1.0;          // P1
                for (int n = 1; n < maxOrder; n++) {              // recurrence
                    P[n + 1] = ((2.0 * n + 1.0) * (2.0 * x - 1.0) * P[n]
                            - n * P[n - 1]) / (n + 1.0);
                }

                // Accumulate the weighted sum: Σ (2n+1) * F[i][n] * Pn
                double sum = 0.0;
                for (int n = 0; n <= maxOrder; n++) {
                    sum += (2.0 * n + 1.0) * F[i][n] * P[n];
                }
                W[i][k] = GW[k] * sum;
            }
        }
    }

    double[] getFractionalOffsets() {
        double[] result = new double[C.length + 1];
        result[0] = 0.0;
        System.arraycopy(C, 0, result, 1, C.length); // Return a copy of the array, for safety
        return result;

    }

}
