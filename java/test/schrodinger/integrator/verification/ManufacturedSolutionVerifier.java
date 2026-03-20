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

    protected ManufacturedSolutionVerifier(DoubleUnaryOperator f, // The exact solution to the problem
                                           DoubleUnaryOperator fPrime, // First derivative of the exact solution
                                           DoubleUnaryOperator qTilde,  // the Q(x) function, - f''(x) / f(x)
                                           DoubleUnaryOperator qTildePrime, // Q'(x)
                                           DoubleUnaryOperator qTildeDoublePrime, // Q''(x)
                                           Integrator integrator) {
        this.f = f;
        this.fPrime = fPrime;
        this.qTilde = qTilde;
        this.qTildePrime = qTildePrime;
        this.qTildeDoublePrime = qTildeDoublePrime;
        this.integrator = integrator;
    }

    public void propagate(double xmin, double xmax, int nPointsStart, int nPointsEnd, int nPointsStep, Integrator.Direction direction, ConvergenceParams params) {
        int d = direction.getValue();

        System.out.println(integrator.getClass().getSimpleName());
        {
            String s = (d == 1) ? "error@xmax" : "error@xmin";
            System.out.printf("%10s %27s %27s %27s", "np", "step", s, "errorRMS\n");
        }

        double maxRatio = 0.;
        for (int np = nPointsStart; np <= nPointsEnd; np += nPointsStep) {
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

            double rmsError = 0.;
            {
                final int startIndex2 = isBackward ? np - integrator.minHistoryLength() : integrator.minHistoryLength() - 1;
                final int endIndex2 = isBackward ? np - 1 - integrator.minHistoryLength() : np - 1;
                Predicate<Integer> condition = isBackward ? n -> (n > 0) : n -> (n < endIndex2);

                double count = 0.;
                for (int n = startIndex2; condition.test(n); n = n + d) {
                    psi[n + d] = integrator.propagate(psi, currentPsiPrime, n, step,
                            qTildeI, qTildePrimeI, qTildeDoublePrimeI,
                            direction);
                    rmsError += Math.pow(psi[n + d] - f.applyAsDouble(xmin + (n + d) * step), 2);
                    count++;
                }

                rmsError = Math.sqrt(rmsError / count);
            }

            double error = isBackward ? f.applyAsDouble(xmin) - psi[0] : f.applyAsDouble(xmax) - psi[np - 1];
            System.out.printf("%10d %27.18f %27.20f %27.20f\n", np, step, error, rmsError);

            String className = integrator.getClass().getSimpleName();
            double errorBound = params.calculateMaxError(step);
            assertTrue(Math.abs(error) < errorBound,
                    String.format("%s : |error| = %.3e should be < %.3e for np = %s",
                            className, Math.abs(error), errorBound, np));
            double ratio = errorBound / Math.abs(error);
            if (ratio > maxRatio) maxRatio = ratio;
        }

        System.out.printf("Max ratio MAX_ERROR / ACTUAL_ERROR = %20.4f\n", maxRatio);
        double maxRatioThreshold = 100;
        if (maxRatio > maxRatioThreshold) {
            System.out.printf("Warning: Ratio is greater than the threshold %20.4f. Consider tightnening "
                    + "the errorBound threshold to strengthen the test!\n", maxRatioThreshold);
        }

    }
}