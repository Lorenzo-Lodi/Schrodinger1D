package schrodinger.integrator.notworking;

import schrodinger.integrator.Integrator;

import java.util.function.IntToDoubleFunction;

public class Burke1980 implements Integrator {

    @Override
    public double propagate(double[] psi, int n, double step, IntToDoubleFunction qTildeFunction, Direction direction) {

        final double h = step;
        final double h2 = h * h;

        int d = direction.getValue();
        int nP3 = n + 3 * d;
        int nP2 = n + 2 * d;
        int nP1 = n + d;
        int n0 = n;
        int n1 = n - d;
        int n2 = n - 2 * d;
        int n3 = n - 3 * d;

        double Q_n2 = qTildeFunction.applyAsDouble(n2);
        double Q_n1 = qTildeFunction.applyAsDouble(n1);
        double Q_n0 = qTildeFunction.applyAsDouble(n0);
        double Q_nP1 = qTildeFunction.applyAsDouble(nP1);
        double Q_nP2 = qTildeFunction.applyAsDouble(nP2);


        return 1. / 465. * (-1920. * psi[n0] + 4770. * psi[n1] - 1920. * psi[n2] - 465. * psi[n3] +
                h2 * (-2538. * Q_n0 - 688. * Q_n1 - 23. * Q_n2 - 688. * Q_nP1 - 23. * Q_nP2));
    }

    @Override
    public int minHistoryLength() {
        return 3;
    }

    @Override
    public int globalConvergenceOrder() {
        return 4;
    }
}
