package schrodinger.integrator.notworking;

import schrodinger.integrator.Integrator;

import java.util.function.IntToDoubleFunction;

public class PadeNumerov8 implements Integrator {

    @Override
    public int minHistoryLength() { return 2; }

    @Override
    public double propagate(double[] psi,
                         int n,
                         double step,
                         IntToDoubleFunction qTildeFunction,
                         Direction direction) {

        int d = direction.getValue();
        double h = step;

        double h2 = h * h;
        double h4 = h2 * h2;
        double h6 = h4 * h2;

        double Qn   = qTildeFunction.applyAsDouble(n);
        double Qnm1 = qTildeFunction.applyAsDouble(n - d);
        double Qnp1 = qTildeFunction.applyAsDouble(n + d);

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
