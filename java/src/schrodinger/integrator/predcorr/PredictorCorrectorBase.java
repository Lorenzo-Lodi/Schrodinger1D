package schrodinger.integrator.predcorr;

import schrodinger.QuantumLevel;
import schrodinger.integrator.Integrator;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialHarmonic;
import schrodinger.potential.SchrodingerSystem;

import java.util.function.IntToDoubleFunction;

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
public abstract class PredictorCorrectorBase implements Integrator {

    protected final Integrator predictor;

    protected PredictorCorrectorBase(Integrator predictor) {
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
     * @param dir   integration direction
     * @return array of length {@code steps} where result[i] = predicted psi at n+(i+1)*d
     */
    protected double[] predictAhead(double[] psi, int n, int steps,
                                    double step, QuantumLevel level, Direction dir) {
        return predictAheadImpl(psi[n], psi, n, steps, step, level, dir);
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
     * @param dir    integration direction
     * @return array of length {@code steps} where result[i] = predicted psi at n+(i+1)*d
     */
    protected double[] predictAhead(double psiAtN, double[] psi, int n, int steps,
                                    double step, QuantumLevel level, Direction dir) {
        return predictAheadImpl(psiAtN, psi, n, steps, step, level, dir);
    }

    private double[] predictAheadImpl(double psiAtN, double[] psi, int n, int steps,
                                      double step, QuantumLevel level, Direction dir) {
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


            // TODO MESSED UP FIX!!!
            final int globalOffset = baseIdx - nib;
            QuantumLevel adjustedLevel = new QuantumLevelOffsetView(level, globalOffset);
            result[i] = predictor.propagate(buf, nib, step, adjustedLevel, dir);
//            final int fb = baseIdx, fn = nib;
//            IntToDoubleFunction adjQFn = idx -> qFn.applyAsDouble(fb - fn + idx);
//            result[i] = predictor.propagate(buf, nib, step, level, dir);
        }
        return result;
    }

    private static class QuantumLevelOffsetView extends QuantumLevel {
        private final QuantumLevel base;
        private final int offset; // globalIdx = localIdx + offset

        QuantumLevelOffsetView(QuantumLevel base, int offset) {
            super();
            this.base = base;
            this.offset = offset;
        }

        // Override ONLY the methods your predictor actually calls.
        // Start with just this one; add others only if you get errors:
        @Override
        public double QTildeAtGridPoint(double localIdx) {
            return base.QTildeAtGridPoint(localIdx + offset);
        }

        // Example for derivatives (add only if predictor uses them):
        // @Override public double dQTildeDrAtGridPoint(int localIdx) {
        //     return base.dQTildeDrAtGridPoint(localIdx + offset);
        // }
    }

}
