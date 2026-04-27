package schrodinger.integrator;

import java.util.function.DoubleUnaryOperator;
import java.util.function.Predicate;

public class StabilityTest {

    public static void main(String[] args) {

        DoubleUnaryOperator f = x -> Math.exp(-x); // The exact solution to the problem
        DoubleUnaryOperator fPrime = x -> -Math.exp(-x); // First derivative of the exact solution
        DoubleUnaryOperator qTilde = x -> -1.;  // the Q(x) function, - f''(x) / f(x)
        DoubleUnaryOperator qTildePrime = x -> 0.; // Q'(x)
        DoubleUnaryOperator qTildeDoublePrime = x -> 0.; // Q''(x)

        Integrator integrator = IntegratorFactory.getNumerov();
        String className = integrator.getClass().getSimpleName();
        System.out.println(className);


        final double xmin = 0.;
        final double xmax = 40.;
        final int np = 100;
        double step = (xmax - xmin) / (np - 1);
        double[] psi = new double[np];
        double[] currentPsiPrime = new double[1];

        DoubleUnaryOperator qTildeI = i -> qTilde.applyAsDouble(xmin + i * step);
        DoubleUnaryOperator qTildePrimeI = i -> qTildePrime.applyAsDouble(xmin + i * step);
        DoubleUnaryOperator qTildeDoublePrimeI = i -> qTildeDoublePrime.applyAsDouble(xmin + i * step);

        // Initialize
        for (int n = 0; n < integrator.minHistoryLength(); n++) {
            double x = xmin + n * step;
            psi[n] = f.applyAsDouble(x);
            currentPsiPrime[0] = fPrime.applyAsDouble(x);
        }


    }


}
