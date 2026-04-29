package schrodinger.integrator.misc;

import schrodinger.integrator.Integrator;
import schrodinger.integrator.pt_correction.PTCorrector;
import schrodinger.integrator.pt_correction.VerletPTCorrector;

import java.util.function.DoubleUnaryOperator;

/**
 * Verlet method for integrating the Schrödinger equation.
 * Second-order accurate method.
 */
public class Verlet implements Integrator {

    @Override
    public int minHistoryLength() {
        return 2;
    }

    @Override
    public int globalConvergenceOrder() {
        return 2;
    }

    @Override
    public PTCorrector getPertubativeCorrector() {
        return new VerletPTCorrector();
    }

    @Override
    public double propagate(double[] psi, double[] currentPsiPrime, int n, double step,
                            DoubleUnaryOperator qTilde,
                            DoubleUnaryOperator qTildePrime,
                            DoubleUnaryOperator qTildeDoublePrime, Direction direction) {
        return psi[n] * (2.0d - step * step * qTilde.applyAsDouble(n)) - psi[n - direction.getValue()];
    }

}
