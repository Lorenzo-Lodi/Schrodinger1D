package schrodinger.integrator.notworking;

import schrodinger.QuantumState;
import schrodinger.integrator.Integrator;

public class Cowell8 implements Integrator {

    @Override
    public double propagate(double[] psi,
                            int n,
                            QuantumState state,
                            Direction direction) {

        int d = direction.getValue();
        double h = state.getGrid().getStepSizeYCoordinate();
        double h2 = h * h;
        double h4 = h2 * h2;

        double Qn   = state.QTildeValueAt(n);
        double Qnm  = state.QTildeValueAt(n - d);
        double Qnp  = state.QTildeValueAt(n + d);

        double yn   = psi[n];
        double ynm  = psi[n - d];

        double numerator =
                2.0 * yn *
                        (1.0
                                - (5.0/6.0) * h2 * Qn
                                + (1.0/24.0) * h4 * Qn * Qn)
                        -
                        ynm *
                                (1.0
                                        + (1.0/12.0) * h2 * Qnm
                                        - (1.0/360.0) * h4 * Qnm * Qnm);

        double denominator =
                1.0
                        + (1.0/12.0) * h2 * Qnp
                        - (1.0/360.0) * h4 * Qnp * Qnp;

        return numerator / denominator;
    }
}
