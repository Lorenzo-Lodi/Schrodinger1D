package schrodinger.integrator;

import java.util.function.IntToDoubleFunction;

/**
 * Predictor-Corrector Störmer method of order 8 for integrating the Schrödinger equation.
 * <p>
 * Uses an injected predictor (default: Numerov) to estimate future psi values
 * (psi[n+1], psi[n+2], psi[n+3]), then applies a symmetric high-order corrector.
 * <p>
 * The key insight: since Q is known analytically everywhere, predicted values
 * only appear as f = -Q·psi_predicted. The prediction error enters one order
 * lower and does NOT degrade the corrector's order (standard PECE result).
 * <p>
 * Available correctors (all using y_{n+1} - 2y_n + y_{n-1} = h²·Σβ_k·f_{n+k}):
 * <p>
 * ORDER 8  — symmetric {-3..+3}:
 * β: ±3→ 31/60480,  ±2→ -73/10080,  ±1→ 2171/20160,  0→ 12067/15120
 * Error constant: -289/3628800 ≈ 8.0e-5  (52× smaller than Numerov)
 * Needs: psi[n-1], psi[n-2], psi[n-3] as history; predicts psi[n+1..n+3]
 */
public class PredictorCorrector8 extends PredictorCorrectorBase {

    private static final double B3 =    31.0 / 60480.0;
    private static final double B2 =   -73.0 / 10080.0;
    private static final double B1 =  2171.0 / 20160.0;
    private static final double B0 = 12067.0 / 15120.0;

    public PredictorCorrector8() { super(new Numerov()); }
    public PredictorCorrector8(Integrator predictor) { super(predictor); }

    @Override
    public int minHistoryLength() { return 4; }

    @Override
    public double propagate(double[] psi, int n, double step, IntToDoubleFunction qTildeFunction, Direction direction) {
        double h = step;
        double h2 = h * h;
        int d = direction.getValue();

        int nP3 = n + 3 * d;
        int nP2 = n + 2 * d;
        int nP1 = n + d;
        int n0  = n;
        int n1  = n - d;
        int n2  = n - 2 * d;
        int n3  = n - 3 * d;

        double Q_n3  = qTildeFunction.applyAsDouble(n3);
        double Q_n2  = qTildeFunction.applyAsDouble(n2);
        double Q_n1  = qTildeFunction.applyAsDouble(n1);
        double Q_n0  = qTildeFunction.applyAsDouble(n0);
        double Q_nP1 = qTildeFunction.applyAsDouble(nP1);
        double Q_nP2 = qTildeFunction.applyAsDouble(nP2);
        double Q_nP3 = qTildeFunction.applyAsDouble(nP3);

        // ── Predict psi[n+1], psi[n+2], psi[n+3] using the injected predictor ──
        double[] pred = predictAhead(psi, n, 3, step, qTildeFunction, direction);
        double psi_nP2_pred = pred[1];
        double psi_nP3_pred = pred[2];

        // ── Correct using symmetric {-3..+3} formula ──
        // β: {±3: 31/60480,  ±2: -73/10080,  ±1: 2171/20160,  0: 12067/15120}
        double rhs = 2.0 * psi[n0]
                - (1.0 + h2 * B1 * Q_n1) * psi[n1]
                - h2 * (B0 * Q_n0 * psi[n0]
                + B2 * Q_nP2 * psi_nP2_pred + B2 * Q_n2 * psi[n2]
                + B3 * Q_nP3 * psi_nP3_pred + B3 * Q_n3 * psi[n3]);

        return rhs / (1.0 + h2 * B1 * Q_nP1);
    }

}
