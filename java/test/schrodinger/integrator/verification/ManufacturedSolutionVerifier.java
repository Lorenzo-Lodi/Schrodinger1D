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
        int[] points = gridParamsToArray(nPointsStart, nPointsEnd, nPointsStep);
        propagate(xmin, xmax, points, direction, params);
    }

    static int[] gridParamsToArray(int nPointsStart, int nPointsEnd, int nPointsStep) {
        int size = 1 + (nPointsEnd - nPointsStart) / nPointsStep;
        int[] points = new int[size];
        for (int i = 0; i < size; i++) {
            points[i] = nPointsStart + nPointsStep * i;
        }
        return points;
    }

    public void propagate(double xmin, double xmax, int[] points, Integrator.Direction direction, ConvergenceParams params) {
        int d = direction.getValue();
        double t0 = System.nanoTime();
        String className = integrator.getClass().getSimpleName();
        System.out.println(className);
        {
            String s = (d == 1) ? "error@xmax" : "error@xmin";
            System.out.printf("%10s %27s %27s %27s", "np", "step", s, "errorRMS\n");
        }

        double maxRatio = 0.;
        for (int np : points) {
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
            String fmtErrorExtreme = (Math.abs(error) < 1000.) ? "%27.20f" : "%27.10e";
            String fmtErrorRMS = (Math.abs(rmsError) < 1000.) ? "%27.20f" : "%27.10e";
            System.out.printf("%10d %27.18f " + fmtErrorExtreme + " " + fmtErrorRMS + "\n", np, step, error, rmsError);

            double errorBound = params.calculateMaxError(step);
            assertTrue(Math.abs(error) < errorBound,
                    String.format("%s : |error| = %.3e should be < %.3e for np = %s",
                            className, Math.abs(error), errorBound, np));
            double ratio = errorBound / Math.abs(error);
            if (ratio > maxRatio) maxRatio = ratio;
        }

        {
            String fmt = (maxRatio < 1e5) ? "%20.4f" : "%20.6e";
            System.out.printf("Max ratio MAX_ERROR / ACTUAL_ERROR = " + fmt + "\n", maxRatio);
        }
        double elapsedMs = (System.nanoTime() - t0) * 1e-6;
        System.out.printf("Elapsed time for %s in ms = %10.3f\n", className, elapsedMs);

        double maxRatioThreshold = 25000;
        if (maxRatio > maxRatioThreshold) {
            System.out.printf("Ratio MAX_ERROR / ACTUAL_ERROR is greater than the maximum threshold %20.4f. " +
                    "This means the test is too permissive. Tighten the errorBound threshold to strengthen the test!\n", maxRatioThreshold);
        }

    }
}