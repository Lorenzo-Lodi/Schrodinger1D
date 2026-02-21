package schrodinger.integrator;

import java.util.function.IntToDoubleFunction;

/**
 * Abstract base class for predictor-corrector integrators.
 * <p>
 * Holds an injected {@link Integrator} used as the predictor (default: Numerov),
 * and provides {@code predictAhead} helpers that build a chain of predicted psi
 * values without modifying the main psi array.
 * <p>
 * Note on predictor compatibility: the 3-element buffer approach used internally
 * is designed for 2-point predictors (those reading only psi[n] and psi[n-1]).
 * Numerov and TaylorThreePoints both qualify. Multi-point predictors like Stormer5
 * would require a wider buffer and are not suitable as predictors here.
 */
public abstract class PredictorCorrectorBase implements Integrator {

    protected final Integrator predictor;

    protected PredictorCorrectorBase(Integrator predictor) {
        if (predictor.minHistoryLength() != 2) {
            throw new IllegalArgumentException(
                predictor.getClass().getSimpleName()
                + " has minHistoryLength() = " + predictor.minHistoryLength()
                + "; only 2-point predictors are compatible with predictAhead()");
        }
        this.predictor = predictor;
    }

    /**
     * Builds a prediction chain of {@code steps} values starting from position
     * {@code n} in the grid, using the injected predictor.
     * Works for both FORWARD and BACKWARD directions.
     *
     * @param psi   the wavefunction array (read-only; psi[n] and psi[n-d] are used as seeds)
     * @param n     starting grid index (psi[n] is the "current" seed)
     * @param steps number of predicted values to produce
     * @param step  grid step size h
     * @param qFn   Q-tilde function over grid indices
     * @param dir   integration direction
     * @return array of length {@code steps} where result[i] = predicted psi at n+(i+1)*d
     */
    protected double[] predictAhead(double[] psi, int n, int steps,
                                    double step, IntToDoubleFunction qFn, Direction dir) {
        return predictAhead(psi[n], psi[n - dir.getValue()], n, steps, step, qFn, dir);
    }

    /**
     * Overload accepting explicit psi_curr / psi_prev starting values,
     * needed when the prediction chain starts from a corrected value rather
     * than directly from the psi array.
     *
     * @param psi_curr value at grid index {@code n}
     * @param psi_prev value at grid index {@code n - d}
     * @param n        starting grid index (where psi_curr lives)
     * @param steps    number of predicted values to produce
     * @param step     grid step size h
     * @param qFn      Q-tilde function over grid indices
     * @param dir      integration direction
     * @return array of length {@code steps} where result[i] = predicted psi at n+(i+1)*d
     */
    protected double[] predictAhead(double psi_curr, double psi_prev,
                                    int n, int steps,
                                    double step, IntToDoubleFunction qFn, Direction dir) {
        int d = dir.getValue();
        double[] result = new double[steps];
        double prev = psi_prev;
        double curr = psi_curr;

        for (int i = 0; i < steps; i++) {
            final int baseIdx = n + i * d;
            final double fPrev = prev, fCurr = curr;
            // Minimal 3-element buffer: [prev, curr, prediction-slot]
            // Layout flips for BACKWARD so the predictor's direction logic still works.
            double[] buf = (d == 1)
                    ? new double[]{fPrev, fCurr, 0.0}
                    : new double[]{0.0,   fCurr, fPrev};
            // Map local buffer indices 0,1,2  →  grid indices baseIdx-1, baseIdx, baseIdx+1
            IntToDoubleFunction adjQFn = idx -> qFn.applyAsDouble(baseIdx - 1 + idx);

            double predicted = predictor.propagate(buf, 1, step, adjQFn, dir);
            result[i] = predicted;
            prev = curr;
            curr = predicted;
        }
        return result;
    }
}
