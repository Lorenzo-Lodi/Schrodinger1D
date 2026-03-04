package schrodinger.integrator.expfitted;

import schrodinger.QuantumLevel;
import schrodinger.integrator.Integrator;

import java.util.function.IntToDoubleFunction;

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
    public double propagate(double[] psi, int n, double step, QuantumLevel level, Direction direction) {
        double h2 = step * step;
        return (psi[n] * (2.0d - 5. * h2 * level.QTildeAtGridPoint(n) / 6.) -
                (1. + h2 * level.QTildeAtGridPoint(n - direction.getValue()) / 12.)
                        * psi[n - direction.getValue()]) /
                (1. + h2 * level.QTildeAtGridPoint(n + direction.getValue()) / 12.);
    }

}
