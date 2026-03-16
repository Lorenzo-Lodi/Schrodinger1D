package schrodinger.integrator.verification;

import schrodinger.integrator.Integrator;

public class ExpMinusXVerifier {

    public static void propagate_forward_exp_minus_x(Integrator integrator) {
        // NB don't use too many points because some methods have fairly bad roundoff error which distorts the convergence patterns
        for (int np = 10; np <= 20; np++) {
            double xmin = 0.;
            double xmax = 1.;
            double[] psi = new double[np];
            double[] currentPsiPrime = new double[1];
            double step = (xmax - xmin) / (np - 1);

            for (int i = 0; i < integrator.minHistoryLength(); i++) {
                double x = xmin + i * step;
                psi[i] = Math.exp(-x);
                currentPsiPrime[0] = -Math.exp(-x);
            }

            for (int n = integrator.minHistoryLength() - 1; n < np - 1; n++) {
                psi[n + 1] = integrator.propagate(psi, currentPsiPrime, n, step,
                        i -> -1., i -> 0., i -> 0.,
                        Integrator.Direction.FORWARD);
            }

            double error = Math.exp(-1.) - psi[np - 1];
            System.out.println(np + " " + error);
        }
    }
}
