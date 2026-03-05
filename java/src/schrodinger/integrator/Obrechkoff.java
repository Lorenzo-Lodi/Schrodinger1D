package schrodinger.integrator;

import java.util.function.IntToDoubleFunction;

public class Obrechkoff implements Integrator {

    @Override
    public int minHistoryLength() {
        return 2;
    }

    @Override
    public int globalConvergenceOrder() {
        return 6;
    }

    @Override
    public double propagate(double[] psi, int n, double step, IntToDoubleFunction qTilde, Direction direction) {
        double eps = 1e-6;
        // TODO !!!!
        return psi[n] * (2.0d - step * step * qTilde.applyAsDouble(n)) - psi[n - direction.getValue()];
    }

}
