package schrodinger.integrator;

import java.util.List;
import java.util.function.DoubleUnaryOperator;

public class StabilityTest {

    public static void main(String[] args) {

//        DoubleUnaryOperator f = x -> Math.exp(-x); // The exact solution to the problem
//        DoubleUnaryOperator fPrime = x -> -Math.exp(-x); // First derivative of the exact solution
//        DoubleUnaryOperator qTilde = x -> -1.;  // the Q(x) function, - f''(x) / f(x)
//        DoubleUnaryOperator qTildePrime = x -> 0.; // Q'(x)
//        DoubleUnaryOperator qTildeDoublePrime = x -> 0.; // Q''(x)
//
        // Let us use a perturbed decaying exponential to avoid special cases.
        double eps = 0.1;
        DoubleUnaryOperator f = x -> Math.exp(-x - eps * Math.sin(x)); // The exact solution to the problem
        DoubleUnaryOperator fPrime = x -> -Math.exp(-x - eps * Math.sin(x)) * (1. + eps * Math.cos(x)); // First derivative of the exact solution
        DoubleUnaryOperator qTilde = x -> -Math.pow(1. + eps * Math.cos(x), 2) - eps * Math.sin(x);  // the Q(x) function, - f''(x) / f(x)
        DoubleUnaryOperator qTildePrime = x -> eps * (2. * Math.sin(x) - Math.cos(x) * (1. - 2. * eps * Math.sin(x))); // Q'(x)
        DoubleUnaryOperator qTildeDoublePrime = x -> eps * (2. * Math.cos(x) + 2. * eps * Math.cos(2. * x) + Math.sin(x)); // Q''(x)


//        List<Integrator> integrators = List.of(IntegratorFactory.getObrechkoff6());
        List<Integrator> integrators = IntegratorFactory.getAll();
        for (Integrator integrator : integrators) {
            String className = integrator.getClass().getSimpleName();
            System.out.println(className);


            final double xmin = 0.;
            final double xmax = 40.;
            for (int k = 0; k <= 40; k++) {
                int np = 50 + 5 * k * k;
                double step = (xmax - xmin) / (np - 1);
                double[] psi = new double[np];
                double[] currentPsiPrime = new double[1];

                DoubleUnaryOperator qTildeI = i -> qTilde.applyAsDouble(xmin + i * step);
                DoubleUnaryOperator qTildePrimeI = i -> qTildePrime.applyAsDouble(xmin + i * step);
                DoubleUnaryOperator qTildeDoublePrimeI = i -> qTildeDoublePrime.applyAsDouble(xmin + i * step);

                // Initialize
                int n;
                double x = 0;
                for (n = 0; n < integrator.minHistoryLength(); n++) {
                    x = xmin + n * step;
                    psi[n] = f.applyAsDouble(x);
                    currentPsiPrime[0] = fPrime.applyAsDouble(x);
                }

                for (n = integrator.minHistoryLength() - 1; n < np - 1; n++) {
                    psi[n + 1] = integrator.propagate(psi, currentPsiPrime, n, step,
                            qTildeI, qTildePrimeI, qTildeDoublePrimeI, Integrator.Direction.FORWARD);

                    x = xmin + (n + 1) * step;
                    double ratio = Math.abs(psi[n + 1]) / f.applyAsDouble(x);
                    if (ratio >= 2.) break;
                }
                System.out.println(-Math.log10(step) + " " + x);
//                System.out.println(step + " " + x);
            }
        }

    }

}
