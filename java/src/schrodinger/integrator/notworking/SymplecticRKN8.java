package schrodinger.integrator.notworking;

import schrodinger.integrator.Integrator;

import java.util.function.IntToDoubleFunction;

/**
 * Symplectic Runge-Kutta-Nyström 8th-order method.
 * Uses Yoshida composition for symplectic integration.
 */
public class SymplecticRKN8 implements Integrator {

    // Yoshida 8th-order composition weights
    private static final double w1 = 0.741670364350612953;
    private static final double w2 = -0.409100825800031594;
    private static final double w3 = 0.190754710296238379;
    private static final double w0 = 1.0 - 2.0 * (w1 + w2 + w3);

    private static final double[] W = {
            w1, w2, w3, w0, w3, w2, w1
    };

    @Override
    public int minHistoryLength() { return 2; }

    @Override
    public double propagate(double[] psi,
                         int n,
                         double step,
                         IntToDoubleFunction qTildeFunction,
                         Direction direction) {

        int d = direction.getValue();
        double h = step;

        // Current and previous values
        double yPrev = psi[n - d];
        double yCurr = psi[n];

        double yNext = yCurr;
        double yBefore = yPrev;

        for (double w : W) {

            double hh = w * h;
            double hh2 = hh * hh;

            // Current Q value
            double q = qTildeFunction.applyAsDouble(n);

            // Verlet kernel in 2-point form
            double yTemp = 2.0 * yNext - yBefore + hh2 * q * yNext;

            yBefore = yNext;
            yNext = yTemp;
        }

        return yNext;
    }
}
