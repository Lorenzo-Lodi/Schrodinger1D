package schrodinger.integrator.predcorr;

import schrodinger.integrator.Integrator;
import schrodinger.integrator.expfitted.Numerov;

import java.util.function.DoubleUnaryOperator;


/**
 * Predictor-Corrector Störmer method of order 6 for integrating the Schrödinger equation.
 * <p>
 * Uses an injected predictor (default: Numerov) to estimate future psi values
 * (psi[n+1], psi[n+2]), then applies a symmetric high-order corrector.
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
 */
public class PredictorCorrector6 extends PredictorCorrectorBase {

    public PredictorCorrector6() {
        super(new Numerov());
    }

    public PredictorCorrector6(Integrator predictor) {
        super(predictor);
    }

    @Override
    public int minHistoryLength() {
        return Math.max(3, predictor.minHistoryLength());
    }

    @Override
    public int globalConvergenceOrder() {
        return 6;
    }

    @Override
    public double propagate(double[] psi, double[] currentPsiPrime, int n, double step,
                            DoubleUnaryOperator qTilde,
                            DoubleUnaryOperator qTildePrime,
                            DoubleUnaryOperator qTildeDoublePrime, Direction direction) {
        double h = step;
        double h2 = h * h;
        int d = direction.getValue();

        int nP2 = n + 2 * d;  // n+2
        int nP1 = n + d;      // n+1  (target)
        int n0 = n;            // n
        int n1 = n - d;        // n-1
        int n2 = n - 2 * d;   // n-2

        double Q_n2 = qTilde.applyAsDouble(n2);
        double Q_n1 = qTilde.applyAsDouble(n1);
        double Q_n0 = qTilde.applyAsDouble(n0);
        double Q_nP1 = qTilde.applyAsDouble(nP1);
        double Q_nP2 = qTilde.applyAsDouble(nP2);

        // ── Predict psi[n+1] and psi[n+2] using the injected predictor ──
        double[] pred = predictAhead(psi, n, 2, step, qTilde, direction);
        double psi_nP2_pred = pred[1];

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
