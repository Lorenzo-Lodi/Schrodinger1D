package schrodinger.integrator;

import java.util.function.IntToDoubleFunction;

/**
 * Numerov's method for integrating the Schrödinger equation.
 * Fourth-order accurate method.
 */
public class Numerov implements Integrator {

    @Override
    public double propagate(double[] psi, int n, double step, IntToDoubleFunction qTildeFunction, Direction direction) {
        double h2 = step * step;
        return (psi[n] * (2.0d - 5. * h2 * qTildeFunction.applyAsDouble(n) / 6.) -
                (1. + h2 * qTildeFunction.applyAsDouble(n - direction.getValue()) / 12.)
                        * psi[n - direction.getValue()]) /
                (1. + h2 * qTildeFunction.applyAsDouble(n + direction.getValue()) / 12.);
    }

}
