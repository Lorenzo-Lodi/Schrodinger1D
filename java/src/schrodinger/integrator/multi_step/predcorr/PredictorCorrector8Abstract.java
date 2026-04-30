package schrodinger.integrator.multi_step.predcorr;

import schrodinger.integrator.Integrator;

import java.util.function.DoubleUnaryOperator;

/**
 * Predictor-Corrector Störmer method of order 8 for integrating Schrödinger equation.
 * <p>
 * Uses an injected predictor (default: Numerov) to estimate future psi values,
 * then applies an 8th-order symmetric corrector with stencil {-3..+3}.
 * <p>
 * Implements a P(EC)^N scheme: one initial prediction followed by {@code CORRECTIONS}
 * correction passes, each preceded by a re-prediction of the future values (except
 * the last correction, which uses the values from the previous re-prediction).
 * <p>
 * Coefficients for ODE form y'' = -Q(x)·y:
 * β: ±3→ 31/60480,  ±2→ -73/10080,  ±1→ 2171/20160,  0→ 12067/15120
 * Error constant: -289/3628800 ≈ -8.0e-5
 */
public abstract class PredictorCorrector8Abstract extends PredictorCorrectorBase {

    private static final double B3 = 31.0 / 60480.0;
    private static final double B2 = -73.0 / 10080.0;
    private static final double B1 = 2171.0 / 20160.0;
    private static final double B0 = 12067.0 / 15120.0;

    private final int maxIterations;

    public PredictorCorrector8Abstract(Integrator predictor, int maxIterations) {
        super(predictor);
        this.maxIterations = maxIterations;
    }

    @Override
    public int minHistoryLength() {
        return Math.max(4, predictor.minHistoryLength());
    }

    @Override
    public int globalConvergenceOrder() {
        return 8;
    }

    @Override
    public double propagate(double[] psi, double[] currentPsiPrime, int n, double step,
                            DoubleUnaryOperator qTilde,
                            DoubleUnaryOperator qTildePrime,
                            DoubleUnaryOperator qTildeDoublePrime, Direction direction) {
        final double h = step;
        final double h2 = h * h;
        final int d = direction.getValue();

        int nP3 = n + 3 * d;
        int nP2 = n + 2 * d;
        int nP1 = n + d;
        int n0 = n;
        int n1 = n - d;
        int n2 = n - 2 * d;
        int n3 = n - 3 * d;

        final double Q_n3 = qTilde.applyAsDouble(n3);
        final double Q_n2 = qTilde.applyAsDouble(n2);
        final double Q_n1 = qTilde.applyAsDouble(n1);
        final double Q_n0 = qTilde.applyAsDouble(n0);
        final double Q_nP1 = qTilde.applyAsDouble(nP1);
        final double Q_nP2 = qTilde.applyAsDouble(nP2);
        final double Q_nP3 = qTilde.applyAsDouble(nP3);

        // ── P: Initial prediction (ψ[n+2] and ψ[n+3] are used; ψ[n+1] is discarded) ──
        double[] pred = predictAhead(psi, n, 3, step, qTilde, direction);
        double psi_nP2 = pred[1];
        double psi_nP3 = pred[2];

        // ── P(EC)^CORRECTIONS loop ──
        double psi_nP1 = 0.0;
        for (int iter = 0; iter < maxIterations; iter++) {
            // C: correct ψ[n+1] using current psi_nP2 / psi_nP3
            double rhs = 2.0 * psi[n0]
                    - (1.0 + h2 * B1 * Q_n1) * psi[n1]
                    - h2 * (B0 * Q_n0 * psi[n0]
                    + B2 * Q_nP2 * psi_nP2 + B2 * Q_n2 * psi[n2]
                    + B3 * Q_nP3 * psi_nP3 + B3 * Q_n3 * psi[n3]);
            psi_nP1 = rhs / (1.0 + h2 * B1 * Q_nP1);

            // E: re-predict ψ[n+2] and ψ[n+3] from corrected ψ[n+1], unless this was the last correction
            if (iter < maxIterations - 1) {
                double[] repred = predictAhead(psi_nP1, psi, nP1, 2, step, qTilde, direction);
                psi_nP2 = repred[0];
                psi_nP3 = repred[1];
            }
        }
        return psi_nP1;
    }

}
