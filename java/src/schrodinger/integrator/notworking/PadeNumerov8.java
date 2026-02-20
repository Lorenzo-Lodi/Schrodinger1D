package schrodinger.integrator.notworking;

import schrodinger.QuantumState;
import schrodinger.integrator.Integrator;

public class PadeNumerov8 {

    public double propagate(double[] psi,
                            int n,
                            QuantumState state,
                            Integrator.Direction direction) {

        int d = direction.getValue();
        double h = state.getGrid().getStepSizeYCoordinate();

        double h2 = h * h;
        double h4 = h2 * h2;
        double h6 = h4 * h2;

        double Qn   = state.QTildeValueAt(n);
        double Qnm1 = state.QTildeValueAt(n - d);
        double Qnp1 = state.QTildeValueAt(n + d);

        double yn   = psi[n];
        double ynm1 = psi[n - d];

        double term_n =
                1.0
                        - (5.0/6.0) * h2 * Qn
                        + (1.0/24.0) * h4 * Qn * Qn
                        - (1.0/720.0) * h6 * Qn * Qn * Qn;

        double term_nm1 =
                1.0
                        + (1.0/12.0) * h2 * Qnm1
                        - (1.0/360.0) * h4 * Qnm1 * Qnm1
                        + (1.0/20160.0) * h6 * Qnm1 * Qnm1 * Qnm1;

        double denominator =
                1.0
                        + (1.0/12.0) * h2 * Qnp1
                        - (1.0/360.0) * h4 * Qnp1 * Qnp1
                        + (1.0/20160.0) * h6 * Qnp1 * Qnp1 * Qnp1;

        double numerator =
                2.0 * yn * term_n
                        - ynm1 * term_nm1;

        return numerator / denominator;
    }
}
