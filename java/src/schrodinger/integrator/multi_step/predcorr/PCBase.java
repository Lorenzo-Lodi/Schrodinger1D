package schrodinger.integrator.multi_step.predcorr;

import schrodinger.integrator.Integrator;

import java.util.function.DoubleUnaryOperator;

/**
 * Abstract base class for predictor-corrector integrators.
 * <p>
 * Holds an injected {@link Integrator} used as the predictor (default: Numerov),
 * and provides {@code predictAhead} helpers that build a chain of predicted psi
 * values without modifying the main psi array.
 * <p>
 * The dynamic sliding-window buffer supports any predictor regardless of its
 * history length, keyed off {@code predictor.minHistoryLength()}.
 */
public abstract class PCBase implements Integrator {

    protected final Integrator predictor;

    protected PCBase(Integrator predictor) {
        this.predictor = predictor;
    }

    /**
     * Builds a prediction chain of {@code steps} values starting from position
     * {@code n} in the grid, using the injected predictor.
     * Works for both FORWARD and BACKWARD directions.
     *
     * @param psi   the wavefunction array (read-only; deep history is read as needed)
     * @param n     starting grid index (psi[n] is the "current" seed)
     * @param steps number of predicted values to produce
     * @param step  grid step size h
     * @param qFn   Q-tilde function over grid indices
     * @param dir   integration direction
     * @return array of length {@code steps} where result[i] = predicted psi at n+(i+1)*d
     */
    protected double[] predictAhead(double[] psi, int n, int steps,
                                    double step, DoubleUnaryOperator qFn, Direction dir) {
        return predictAheadImpl(psi[n], psi, n, steps, step, qFn, dir);
    }

    /**
     * Overload accepting an explicit corrected value for grid index {@code n},
     * needed when the prediction chain starts from a corrected value rather
     * than directly from the psi array.
     *
     * @param psiAtN corrected value at grid index {@code n} (overrides psi[n])
     * @param psi    the wavefunction array (read-only; deep history is read as needed)
     * @param n      starting grid index (where psiAtN lives)
     * @param steps  number of predicted values to produce
     * @param step   grid step size h
     * @param qFn    Q-tilde function over grid indices
     * @param dir    integration direction
     * @return array of length {@code steps} where result[i] = predicted psi at n+(i+1)*d
     */
    protected double[] predictAhead(double psiAtN, double[] psi,
                                    int n, int steps,
                                    double step, DoubleUnaryOperator qFn, Direction dir) {
        return predictAheadImpl(psiAtN, psi, n, steps, step, qFn, dir);
    }

    private double[] predictAheadImpl(double psiAtN, double[] psi, int n, int steps,
                                      double step, DoubleUnaryOperator qFn, Direction dir) {
        int d = dir.getValue();
        int k = predictor.minHistoryLength();
        int nib = (d == 1) ? k - 1 : 1;        // n_in_buf
        double[] buf = new double[k + 1];
        double[] result = new double[steps];

        for (int i = 0; i < steps; i++) {
            int baseIdx = n + i * d;
            for (int j = 0; j <= k; j++) {
                int g = baseIdx - nib + j;
                if (g == baseIdx + d) {
                    buf[j] = 0.0;
                    continue;
                }   // future slot
                if (g == n) {
                    buf[j] = psiAtN;
                    continue;
                } // override
                buf[j] = (d == 1)
                        ? ((g <= n) ? psi[g] : result[g - n - 1])
                        : ((g >= n) ? psi[g] : result[n - g - 1]);
            }
            final int fb = baseIdx, fn = nib;
            DoubleUnaryOperator adjQFn = idx -> qFn.applyAsDouble(fb - fn + idx);
            result[i] = predictor.propagate(buf, null, nib, step, adjQFn, null, null, dir);
        }
        return result;
    }

    @Override
    public boolean needsPotentialCapping() {
        return true;
    }
}
