package schrodinger.integrator.verification;

import schrodinger.integrator.Integrator;

import java.util.function.DoubleUnaryOperator;
import java.util.function.Predicate;

import static org.junit.jupiter.api.Assertions.assertTrue;

public abstract class ManufacturedSolutionVerifier {
    private final DoubleUnaryOperator f;
    private final DoubleUnaryOperator fPrime;
    private final DoubleUnaryOperator qTilde;
    private final DoubleUnaryOperator qTildePrime;
    private final DoubleUnaryOperator qTildeDoublePrime;
    private final Integrator integrator;

    protected ManufacturedSolutionVerifier(DoubleUnaryOperator f,
                                           DoubleUnaryOperator fPrime,
                                           DoubleUnaryOperator qTilde,
                                           DoubleUnaryOperator qTildePrime,
                                           DoubleUnaryOperator qTildeDoublePrime,
                                           Integrator integrator) {
        this.f = f;
        this.fPrime = fPrime;
        this.qTilde = qTilde;
        this.qTildePrime = qTildePrime;
        this.qTildeDoublePrime = qTildeDoublePrime;
        this.integrator = integrator;
    }

    public void propagate(Integrator.Direction direction, ConvergenceParams params) {
        System.out.println(integrator.getClass().getSimpleName());
        int d = direction.getValue();

        for (int np = 100; np <= 500; np += 10) {
            double xmin = 0.;
            double xmax = 4. * Math.PI;
            double[] psi = new double[np];
            double[] currentPsiPrime = new double[1];
            double step = (xmax - xmin) / (np - 1);

            boolean isBackward = direction.equals(Integrator.Direction.BACKWARD);
            {
                final int startIndex1 = isBackward ? np - 1 : 0;
                final int endIndex1 = isBackward ? np - 1 - integrator.minHistoryLength() : integrator.minHistoryLength();
                Predicate<Integer> condition = isBackward ? n -> (n > endIndex1) : n -> (n < endIndex1);
                for (int n = startIndex1; condition.test(n); n = n + d) {
                    double x = xmin + n * step;
                    psi[n] = f.applyAsDouble(x);
                    currentPsiPrime[0] = fPrime.applyAsDouble(x);
                }
            }

            DoubleUnaryOperator qTildeI = i -> qTilde.applyAsDouble(xmin + i * step);
            DoubleUnaryOperator qTildePrimeI = i -> qTildePrime.applyAsDouble(xmin + i * step);
            DoubleUnaryOperator qTildeDoublePrimeI = i -> qTildeDoublePrime.applyAsDouble(xmin + i * step);

            {
                final int startIndex2 = isBackward ? np - integrator.minHistoryLength() : integrator.minHistoryLength() - 1;
                final int endIndex2 = isBackward ? np - 1 - integrator.minHistoryLength() : np - 1;
                Predicate<Integer> condition = isBackward ? n -> (n > 0) : n -> (n < endIndex2);

                for (int n = startIndex2; condition.test(n); n = n + d) {
                    psi[n + d] = integrator.propagate(psi, currentPsiPrime, n, step,
                            qTildeI, qTildePrimeI, qTildeDoublePrimeI,
                            direction);
                }
            }

            double error = isBackward ? f.applyAsDouble(xmin) - psi[0] : f.applyAsDouble(xmax) - psi[np - 1];
            System.out.println(np + " " + error);

            String className = integrator.getClass().getSimpleName();
            double errorBound = params.calculateMaxError(np);
            assertTrue(Math.abs(error) < errorBound,
                    String.format("%s : |error| = %.3e should be < %.3e for np = %s",
                            className, Math.abs(error), errorBound, np));
        }

    }
}