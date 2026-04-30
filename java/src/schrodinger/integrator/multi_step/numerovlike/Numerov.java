package schrodinger.integrator.multi_step.numerovlike;

import schrodinger.integrator.Integrator;
import schrodinger.integrator.pt_correction.NumerovPTCorrector;
import schrodinger.integrator.pt_correction.PTCorrector;

import java.util.function.DoubleUnaryOperator;

/**
 * Numerov's method for integrating the Schrödinger equation.
 * Fourth-order accurate method.
 */
public class Numerov implements Integrator {

    @Override
    public int minHistoryLength() {
        return 2;
    }

    @Override
    public int globalConvergenceOrder() {
        return 4;
    }

    @Override
    public PTCorrector getPertubativeCorrector() {
        return new NumerovPTCorrector();
    }

    @Override
    public boolean needsPotentialCapping() {
        return true;
    }

    @Override
    public double propagate(double[] psi, double[] currentPsiPrime, int n, double step,
                            DoubleUnaryOperator qTilde,
                            DoubleUnaryOperator qTildePrime,
                            DoubleUnaryOperator qTildeDoublePrime, Direction direction) {
        double h2 = step * step;
        return (psi[n] * (2.0d - 5. * h2 * qTilde.applyAsDouble(n) / 6.) -
                (1. + h2 * qTilde.applyAsDouble(n - direction.getValue()) / 12.)
                        * psi[n - direction.getValue()]) /
                (1. + h2 * qTilde.applyAsDouble(n + direction.getValue()) / 12.);
    }

}
