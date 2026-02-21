package schrodinger.integrator;

import java.util.function.IntToDoubleFunction;

/**
 * Predictor-Corrector Störmer method of order 8 for integrating Schrödinger equation.
 * <p>
 * Uses an injected predictor (default: Numerov) to estimate future psi values,
 * then applies an 8th-order symmetric corrector with stencil {-3..+3}.
 * <p>
 * Implements full PECE scheme to achieve 8th-order convergence:
 * Predict -> Evaluate -> Correct -> Evaluate -> Correct
 * <p>
 * Coefficients for ODE form y'' = -Q(x)·y:
 *   β: ±3→ 31/60480,  ±2→ -73/10080,  ±1→ 2171/20160,  0→ 12067/15120
 *   Error constant: -289/3628800 ≈ -8.0e-5
 */
public class PredictorCorrector8Alt extends PredictorCorrectorBase {

    private static final double B3 =    31.0 / 60480.0;
    private static final double B2 =   -73.0 / 10080.0;
    private static final double B1 =  2171.0 / 20160.0;
    private static final double B0 = 12067.0 / 15120.0;

    public PredictorCorrector8Alt() { super(new Numerov()); }
    public PredictorCorrector8Alt(Integrator predictor) { super(predictor); }

    @Override
    public int minHistoryLength() { return 4; }

    @Override
    public double propagate(double[] psi, int n, double step, IntToDoubleFunction qTildeFunction, Direction direction) {
        final double h  = step;
        final double h2 = h * h;
        final int d = direction.getValue();

        final int nP3 = n + 3*d, nP2 = n + 2*d, nP1 = n + d;
        final int n0  = n;
        final int n1  = n - d, n2 = n - 2*d, n3 = n - 3*d;

        final double Q_n3  = qTildeFunction.applyAsDouble(n3);
        final double Q_n2  = qTildeFunction.applyAsDouble(n2);
        final double Q_n1  = qTildeFunction.applyAsDouble(n1);
        final double Q_n0  = qTildeFunction.applyAsDouble(n0);
        final double Q_nP1 = qTildeFunction.applyAsDouble(nP1);
        final double Q_nP2 = qTildeFunction.applyAsDouble(nP2);
        final double Q_nP3 = qTildeFunction.applyAsDouble(nP3);

        // ── P: First prediction pass (from n) ──
        double[] pred = predictAhead(psi, n, 3, step, qTildeFunction, direction);
        double psi_nP2_pred = pred[1];
        double psi_nP3_pred = pred[2];

        // ── C: Correct ψ[n+1] ──
        double rhs = 2.0 * psi[n0]
                - (1.0 + h2 * B1 * Q_n1) * psi[n1]
                - h2 * ( B0 * Q_n0  * psi[n0]
                + B2 * Q_nP2 * psi_nP2_pred + B2 * Q_n2 * psi[n2]
                + B3 * Q_nP3 * psi_nP3_pred + B3 * Q_n3 * psi[n3] );

        double psi_nP1 = rhs / (1.0 + h2 * B1 * Q_nP1);

        // ── E & C: Re-predict ψ[n+2], ψ[n+3] using corrected ψ[n+1], then correct ──
        // Second prediction pass starts from the corrected psi_nP1 as its seed.
        double[] repred = predictAhead(psi_nP1, psi[n0], nP1, 2, step, qTildeFunction, direction);
        double psi_nP2_repred = repred[0];
        double psi_nP3_repred = repred[1];

        // Re-correct ψ[n+1] with improved future values
        double rhs2 = 2.0 * psi[n0]
                - (1.0 + h2 * B1 * Q_n1) * psi[n1]
                - h2 * ( B0 * Q_n0  * psi[n0]
                + B2 * Q_nP2 * psi_nP2_repred + B2 * Q_n2 * psi[n2]
                + B3 * Q_nP3 * psi_nP3_repred + B3 * Q_n3 * psi[n3] );

        return rhs2 / (1.0 + h2 * B1 * Q_nP1);
    }

}
