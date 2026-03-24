package schrodinger.integrator.rungekutta;

import schrodinger.integrator.Integrator;

import java.util.function.DoubleUnaryOperator;

public abstract class CFMagnusAbstract {

    private final int NODES;
    private final int EXPONENTIALS;

    private final double[] C; // Gauss-Legendre nodes
    private final double[][] W; // Weights
    private final double[] A1W;


    public CFMagnusAbstract(int nodes, int exponentials, double[] c, double[][] w, double[] a1W) {
        NODES = nodes;
        EXPONENTIALS = exponentials;
        C = c;
        W = w;
        A1W = a1W;
    }

    public double propagate(double[] psi, double[] currentPsiPrime, int n, double step,
                            DoubleUnaryOperator qTilde,
                            DoubleUnaryOperator qTildePrime,
                            DoubleUnaryOperator qTildeDoublePrime,
                            Integrator.Direction direction) {

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


}
