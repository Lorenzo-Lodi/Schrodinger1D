package schrodinger.integrator;

import org.junit.jupiter.api.Test;

import java.util.function.DoubleUnaryOperator;
import java.util.function.Predicate;

import static org.junit.jupiter.api.Assertions.assertTrue;

public abstract class AbstractIntegratorTest {

    /**
     * Returns the integrator to be tested.
     * Subclasses must implement this method to provide their specific integrator.
     */
    protected abstract Integrator getIntegrator();

    @Test
    void harmonic_oscillator_ground_state() {
        HarmonicOscillatorConvergenceVerifier ho = new HarmonicOscillatorConvergenceVerifier(getIntegrator());
        ho.harmonic_oscillator_ground_state();
    }

    @Test
    void harmonic_oscillator_10th_excited_state() {
        HarmonicOscillatorConvergenceVerifier ho = new HarmonicOscillatorConvergenceVerifier(getIntegrator());
        ho.harmonic_oscillator_10th_excited_state();
    }


    // This is a freeze/regression test to detect accidental changes to integrator behavior during refactoring.
    // Uses "meaningless" values: cos(k) for psi and a polynomial for Q.
    // Does NOT guarantee correctness - other tests verify that.
    // Subclasses implement getIntegrator() to provide specific integrator.
    private double integrateOneStep(double qFactor, Integrator.Direction direction) {
        int arraySize = 20;
        double[] psi = new double[arraySize];
        double[] currentPsiPrime = new double[1];
        for (int k = 0; k < arraySize; k++) {
            psi[k] = Math.cos(k);
            currentPsiPrime[0] = -Math.sin(k);
        }

        Integrator integrator = getIntegrator();
        DoubleUnaryOperator q = i -> qFactor * (0.1111 + 0.2222 * i + 0.3333 * i * i + 0.44444 * i * i * i);
        DoubleUnaryOperator q1 = i -> qFactor * (0.2222 + 2 * 0.3333 * i + 3. * 0.44444 * i * i);
        DoubleUnaryOperator q2 = i -> qFactor * (2 * 0.3333 + 3. * 2. * 0.44444 * i);
        return integrator.propagate(psi, currentPsiPrime, arraySize / 2, 0.123, q, q1, q2, direction);
    }

    public double integrateOneStep(Integrator.Direction direction) {
        return integrateOneStep(1., direction);
    }

    // Same as integrateOneStep but with negated Q to test Z < 0 branch in exponentially fitted methods.
    public double integrateOneStepNegativeQ(Integrator.Direction direction) {
        return integrateOneStep(-1., direction);
    }

    // Same as integrateOneStep but with very small Q  to test series expansion branch (|Z| < threshold).
    public double integrateOneStepSmallZ(Integrator.Direction direction) {
        return integrateOneStep(2.e-5, direction);
    }

    @Test
    public void time_benchmark() {
        IntegratorSpeedBenchmark.time_benchmark(getIntegrator());
    }

    @Test
    public void propagate_forward_exp_minus_x() {
        Integrator integrator = getIntegrator();

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

    @Test
    public void propagate_forward_exp_to_cos_x() {
        propagate_exp_to_cos_x(Integrator.Direction.FORWARD);
    }

    @Test
    public void propagate_backward_exp_to_cos_x() {
        propagate_exp_to_cos_x(Integrator.Direction.BACKWARD);
    }

    public void propagate_exp_to_cos_x(Integrator.Direction direction) {
        Integrator integrator = getIntegrator();
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
                    psi[n] = Math.exp(Math.cos(x));
                    currentPsiPrime[0] = -Math.sin(x) * Math.exp(Math.cos(x));
                }
            }

            DoubleUnaryOperator qTilde = i -> Math.cos(xmin + i * step) - Math.pow(Math.sin(xmin + i * step), 2);
            DoubleUnaryOperator qTildePrime = i -> -Math.sin(xmin + i * step) - 2. * Math.cos(xmin + i * step) * Math.sin(xmin + i * step);
            DoubleUnaryOperator qTildeDoublePrime = i -> -Math.cos(xmin + i * step)
                    - 2. * Math.pow(Math.cos(xmin + i * step), 2)
                    + 2. * Math.pow(Math.sin(xmin + i * step), 2);

            {
                final int startIndex2 = isBackward ? np - integrator.minHistoryLength() : integrator.minHistoryLength() - 1;
                final int endIndex2 = isBackward ? np - 1 - integrator.minHistoryLength() : np - 1;
                Predicate<Integer> condition = isBackward ? n -> (n > 0) : n -> (n < endIndex2);

                for (int n = startIndex2; condition.test(n); n = n + d) {
                    psi[n + d] = integrator.propagate(psi, currentPsiPrime, n, step,
                            qTilde, qTildePrime, qTildeDoublePrime,
                            direction);
                }
            }

            double error = isBackward ? Math.exp(Math.cos(xmin)) - psi[0] : Math.exp(Math.cos(xmax)) - psi[np - 1];
            System.out.println(np + " " + error);

            String className = integrator.getClass().getSimpleName();
            double maxError = error_bound_propagate_exp_to_cos_x(np);
            assertTrue(Math.abs(error) < maxError,
                    String.format("%s : |error| = %.3e should be < %.3e for np = %s",
                            className, Math.abs(error), maxError, np));
        }

    }

    // We assume errors are in the format |error| < (C/np)^a
    // The higher a , the better;
    // The smallest C, the better;
    abstract protected double error_bound_propagate_exp_to_cos_x(int np);


}