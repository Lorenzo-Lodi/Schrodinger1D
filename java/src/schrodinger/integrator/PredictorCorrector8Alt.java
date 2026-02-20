package schrodinger.integrator;

import java.util.function.IntToDoubleFunction;

/**
 * Predictor-Corrector Störmer method of order 8 for integrating Schrödinger equation.
 * <p>
 * Uses Numerov (4th order) as predictor to estimate future psi values,
 * then applies an 8th-order symmetric corrector with stencil {-3..+3}.
 * <p>
 * Implements full PECE scheme to achieve 8th-order convergence:
 * Predict -> Evaluate -> Correct -> Evaluate -> Correct
 * <p>
 * Coefficients for ODE form y'' = -Q(x)·y:
 *   β: ±3→ 31/60480,  ±2→ -73/10080,  ±1→ 2171/20160,  0→ 12067/15120
 *   Error constant: -289/3628800 ≈ -8.0e-5
 */
public class PredictorCorrector8Alt implements Integrator {

    private static final double B3 =    31.0 / 60480.0;
    private static final double B2 =   -73.0 / 10080.0;
    private static final double B1 =  2171.0 / 20160.0;
    private static final double B0 = 12067.0 / 15120.0;

    private double numerovStep(double psi_curr, double psi_prev,
                               double Q_next, double Q_curr, double Q_prev,
                               double h2) {
        double num = (2.0 - 10.0/12.0 * h2 * Q_curr) * psi_curr
                - (1.0 +  1.0/12.0 * h2 * Q_prev) * psi_prev;
        return num / (1.0 + 1.0/12.0 * h2 * Q_next);
    }

    /**
     * Computes the corrector RHS excluding the implicit β₁ term.
     *
     * The corrector formula is:
     *   ψ_{n+1} - 2ψ_n + ψ_{n-1} = h² Σ β_k f_{n+k}
     *
     * where f_{n+k} = -Q_{n+k} ψ_{n+k}.
     *
     * Moving the implicit term (k=+1) to LHS:
     *   ψ_{n+1}(1 + h²β₁Q_{n+1}) = 2ψ_n - ψ_{n-1} + h² Σ_{k≠+1} β_k f_{n+k}
     *
     * This method computes the RHS sum.
     */

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

        // ── P: Predict using Numerov ──
        double psi_nP1_pred = numerovStep(psi[n0], psi[n1], Q_nP1, Q_n0, Q_n1, h2);
        double psi_nP2_pred = numerovStep(psi_nP1_pred, psi[n0], Q_nP2, Q_nP1, Q_n0, h2);
        double psi_nP3_pred = numerovStep(psi_nP2_pred, psi_nP1_pred, Q_nP3, Q_nP2, Q_nP1, h2);

        // ── C: Correct ψ[n+1] ──
        double rhs = 2.0 * psi[n0]
                - (1.0 + h2 * B1 * Q_n1) * psi[n1]
                - h2 * ( B0 * Q_n0  * psi[n0]
                + B2 * Q_nP2 * psi_nP2_pred + B2 * Q_n2 * psi[n2]
                + B3 * Q_nP3 * psi_nP3_pred + B3 * Q_n3 * psi[n3] );

        double psi_nP1 = rhs / (1.0 + h2 * B1 * Q_nP1);

        // ── E & C: Re-predict ψ[n+2] using corrected ψ[n+1], then correct ──
        // Use corrected ψ[n+1] as better "history" for next predictions
        double psi_nP2_repred = numerovStep(psi_nP1, psi[n0], Q_nP2, Q_nP1, Q_n0, h2);
        double psi_nP3_repred = numerovStep(psi_nP2_repred, psi_nP1, Q_nP3, Q_nP2, Q_nP1, h2);

        // Re-correct ψ[n+1] with improved future values
        double rhs2 = 2.0 * psi[n0]
                - (1.0 + h2 * B1 * Q_n1) * psi[n1]
                - h2 * ( B0 * Q_n0  * psi[n0]
                + B2 * Q_nP2 * psi_nP2_repred + B2 * Q_n2 * psi[n2]
                + B3 * Q_nP3 * psi_nP3_repred + B3 * Q_n3 * psi[n3] );

        return rhs2 / (1.0 + h2 * B1 * Q_nP1);
    }

}
