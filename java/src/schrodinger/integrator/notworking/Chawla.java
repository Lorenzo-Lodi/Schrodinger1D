package schrodinger.integrator.notworking;

import schrodinger.integrator.Integrator;

import java.util.function.IntToDoubleFunction;

/**
 * Propagates using a 4-step, 6th-order Extended Numerov method.
 *
 * Requires history: psi[n], psi[n-1], psi[n-2].
 * This method has a local truncation error of O(h^8), making it significantly
 * more accurate than the standard Numerov method (O(h^6)).
 */
public class Chawla implements Integrator {

    @Override
    public int minHistoryLength() { return 3; }

    @Override
    public double propagate(double[] psi, int n, double step, IntToDoubleFunction qTildeFunction, Direction direction) {
        double h = step;
        double h2 = h * h;

        int n_curr  = n;
        int n_prev  = n - direction.getValue();
        int n_prev2 = n - 2 * direction.getValue();
        int n_next  = n + direction.getValue();

        double Q_next  = qTildeFunction.applyAsDouble(n_next);
        double Q_curr  = qTildeFunction.applyAsDouble(n_curr);
        double Q_prev  = qTildeFunction.applyAsDouble(n_prev);
        double Q_prev2 = qTildeFunction.applyAsDouble(n_prev2);

        // Raptis & Cash (1987) asymmetric 6th-order coefficients
        // Stencil offsets: +1, 0, -1, -2
        double b0 =  19.0 / 240.0;   // f_{n+1}
        double b1 =  29.0 /  20.0;   // f_{n}     -- was incorrectly 229/120
        double b2 =   1.0 /   5.0;   // f_{n-1}
        double b3 =  -1.0 / 240.0;   // f_{n-2}   -- was incorrectly +1/120

        // (1 + h²·b0·Q_next)·y_next = 2y_curr - y_prev - h²·(b1·Q_curr·y_curr + b2·Q_prev·y_prev + b3·Q_prev2·y_prev2)
        double term_RHS =
                2.0 * psi[n_curr]
                        - psi[n_prev]
                        - h2 * (  b1 * Q_curr  * psi[n_curr]
                        + b2 * Q_prev  * psi[n_prev]
                        + b3 * Q_prev2 * psi[n_prev2] );

        double denominator = 1.0 + h2 * b0 * Q_next;
        return term_RHS / denominator;
    }
}
