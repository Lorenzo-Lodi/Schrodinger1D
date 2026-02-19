package schrodinger.integrator;

import schrodinger.QuantumState;

public class Stormer7 implements Integrator {

    /**
     * Extended Störmer method for y'' = -Q(x)·y.
     * <p>
     * Uses the recurrence:  y_{n+1} - 2y_n + y_{n-1} = h² · Σ βₖ · f_{n+k}
     * with f = -Q·y (zero-stable, rho(z) = (z-1)^2).
     * <p>
     * Option A (order 5 / LTE O(h^7)): requires psi[n], psi[n-1], psi[n-2], psi[n-3]
     * Option B (order 6 / LTE O(h^8)): requires psi[n], psi[n-1], psi[n-2], psi[n-3], psi[n-4]
     * <p>
     * Coefficients derived from Taylor order conditions; verified symbolically.
     */
    @Override
    public double propagate(double[] psi, int n, QuantumState state, Direction direction) {
        double h = state.getGrid().getStepSizeYCoordinate();
        double h2 = h * h;
        int d = direction.getValue();

        int nP = n + d;        // n+1   (implicit)
        int n0 = n;            // current
        int n1 = n - d;        // n-1
        int n2 = n - 2 * d;    // n-2
        int n3 = n - 3 * d;    // n-3


//         ── Option A: order 5, LTE = O(h^7) ──────────────────────────────────────
//         β: {+1: 19/240,  0: 17/20,  -1: 7/120,  -2: 1/60,  -3: -1/240}

        double b_next = 19.0 / 240.0;
        double b_curr = 17.0 / 20.0;
        double b_prev = 7.0 / 120.0;
        double b_prev2 = 1.0 / 60.0;
        double b_prev3 = -1.0 / 240.0;

        double Q_next = state.QTildeValueAt(nP);
        double Q_curr = state.QTildeValueAt(n0);
        double Q_prev = state.QTildeValueAt(n1);
        double Q_prev2 = state.QTildeValueAt(n2);
        double Q_prev3 = state.QTildeValueAt(n3);

        double rhs = 2.0 * psi[n0] - psi[n1]
                - h2 * (b_curr * Q_curr * psi[n0]
                + b_prev * Q_prev * psi[n1]
                + b_prev2 * Q_prev2 * psi[n2]
                + b_prev3 * Q_prev3 * psi[n3]);
        return rhs / (1.0 + h2 * b_next * Q_next);

    }

}
