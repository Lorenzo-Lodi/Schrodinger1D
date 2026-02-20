package schrodinger.integrator;

import java.util.function.IntToDoubleFunction;

/**
 * Three-point Taylor expansion method for integrating the Schrödinger equation.
 * Second-order accurate method.
 */
public class TaylorThreePoints implements Integrator {

    @Override
    public double propagate(double[] psi, int n, double step, IntToDoubleFunction qTildeFunction, Direction direction) {
        return psi[n] * (2.0d - step * step * qTildeFunction.applyAsDouble(n)) - psi[n - direction.getValue()];
    }

}
