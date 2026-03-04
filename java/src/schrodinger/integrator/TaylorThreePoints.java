package schrodinger.integrator;

import java.util.function.DoubleUnaryOperator;

/**
 * Three-point Taylor expansion method for integrating the Schrödinger equation.
 * Second-order accurate method.
 */
public class TaylorThreePoints implements Integrator {

    @Override
    public int minHistoryLength() { return 2; }

    @Override
    public int globalConvergenceOrder() { return 2; }

    @Override
    public double propagate(double[] psi, int n, double step, DoubleUnaryOperator qTilde, Direction direction) {
        return psi[n] * (2.0d - step * step * qTilde.applyAsDouble(n)) - psi[n - direction.getValue()];
    }

}
