package schrodinger.integrator;

import schrodinger.QuantumState;

public class PredictorCorrector6 implements Integrator {

    /**
     * Predictor-Corrector Störmer methods for y'' = -Q(x)·y.
     * <p>
     * Uses Numerov as predictor to estimate future psi values (psi[n+1], psi[n+2], ...),
     * then applies a symmetric high-order corrector.
     * <p>
     * The key insight: since Q is known analytically everywhere, predicted values
     * only appear as f = -Q·psi_predicted. The prediction error enters one order
     * lower and does NOT degrade the corrector's order (standard PECE result).
     * <p>
     * Available correctors (all using y_{n+1} - 2y_n + y_{n-1} = h²·Σβ_k·f_{n+k}):
     * <p>
     * ORDER 6  — symmetric {-2..+2}:
     * β: ±2→ -1/240,  ±1→ 1/10,  0→ 97/120
     * Error constant: 31/60480 ≈ 5.1e-4   (8× smaller than Numerov)
     * Needs: psi[n-1], psi[n-2] as history; predicts psi[n+1], psi[n+2]
     * <p>
     * ORDER 8  — symmetric {-3..+3}:
     * β: ±3→ 31/60480,  ±2→ -73/10080,  ±1→ 2171/20160,  0→ 12067/15120
     * Error constant: -289/3628800 ≈ 8.0e-5  (52× smaller than Numerov)
     * Needs: psi[n-1], psi[n-2], psi[n-3] as history; predicts psi[n+1..n+3]
     */

// ── Numerov predictor (reusable helper) ──────────────────────────────────────
    private double numerovStep(double psi_curr, double psi_prev,
                               double Q_next, double Q_curr, double Q_prev,
                               double h2) {
        double num = (2.0 - 10.0 / 12.0 * h2 * Q_curr) * psi_curr
                - (1.0 + 1.0 / 12.0 * h2 * Q_prev) * psi_prev;
        return num / (1.0 + 1.0 / 12.0 * h2 * Q_next);
    }

// ── Order 6 predictor-corrector ───────────────────────────────────────────────

    @Override
    public double propagate(double[] psi, int n, QuantumState state, Direction direction) {
        double h = state.getGrid().getStepSizeYCoordinate();
        double h2 = h * h;
        int d = direction.getValue();

        int nP2 = n + 2 * d;  // n+2
        int nP1 = n + d;  // n+1  (target)
        int n0 = n;         // n
        int n1 = n - d;  // n-1
        int n2 = n - 2 * d;  // n-2

        double Q_n2 = state.QTildeValueAt(n2);
        double Q_n1 = state.QTildeValueAt(n1);
        double Q_n0 = state.QTildeValueAt(n0);
        double Q_nP1 = state.QTildeValueAt(nP1);
        double Q_nP2 = state.QTildeValueAt(nP2);

        // ── Predict psi[n+1] and psi[n+2] using Numerov ──
        double psi_nP1_pred = numerovStep(psi[n0], psi[n1], Q_nP1, Q_n0, Q_n1, h2);
        double psi_nP2_pred = numerovStep(psi_nP1_pred, psi[n0], Q_nP2, Q_nP1, Q_n0, h2);

        // ── Correct using symmetric {-2..+2} formula ──
        // β: {±2: -1/240,  ±1: 1/10,  0: 97/120}
        // (1 + h²·(1/10)·Q_{n+1}) · psi[n+1] =
        //     2·psi[n] - (1 + h²·(1/10)·Q_{n-1}) · psi[n-1]
        //   - h²·( (97/120)·Q_n·psi[n]
        //        + (-1/240)·Q_{n+2}·psi[n+2]_pred
        //        + (-1/240)·Q_{n-2}·psi[n-2] )
        double rhs = 2.0 * psi[n0]
                - (1.0 + h2 * (1.0 / 10.0) * Q_n1) * psi[n1]
                - h2 * ((97.0 / 120.0) * Q_n0 * psi[n0]
                + (-1.0 / 240.0) * Q_nP2 * psi_nP2_pred
                + (-1.0 / 240.0) * Q_n2 * psi[n2]);

        return rhs / (1.0 + h2 * (1.0 / 10.0) * Q_nP1);
    }

}
