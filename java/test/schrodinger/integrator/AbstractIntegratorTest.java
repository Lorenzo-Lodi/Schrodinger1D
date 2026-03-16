package schrodinger.integrator;

import org.junit.jupiter.api.Test;

import java.util.function.DoubleUnaryOperator;

public abstract class AbstractIntegratorTest {

    /**
     * Returns the integrator to be tested.
     * Subclasses must implement this method to provide their specific integrator.
     */
    protected abstract Integrator getIntegrator();

    @Test
    void harmonic_oscillator_ground_state() {
        HarmonicOscillatorConvergenceVerifier ho = new HarmonicOscillatorConvergenceVerifier(getIntegrator());
        ho.verify_harmonic_oscillator_ground_state();
    }

    @Test
    void harmonic_oscillator_10th_excited_state() {
        HarmonicOscillatorConvergenceVerifier ho = new HarmonicOscillatorConvergenceVerifier(getIntegrator());
        ho.verify_harmonic_oscillator_10th_excited_state();
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

    protected double integrateOneStep(Integrator.Direction direction) {
        return integrateOneStep(1., direction);
    }

    // Same as integrateOneStep but with negated Q to test Z < 0 branch in exponentially fitted methods.
    protected double integrateOneStepNegativeQ(Integrator.Direction direction) {
        return integrateOneStep(-1., direction);
    }

    // Same as integrateOneStep but with very small Q  to test series expansion branch (|Z| < threshold).
    protected double integrateOneStepSmallZ(Integrator.Direction direction) {
        return integrateOneStep(2.e-5, direction);
    }

    @Test
    void time_benchmark() {
        IntegratorSpeedBenchmark.time_benchmark(getIntegrator());
    }

    @Test
    void propagate_forward_exp_minus_x() {
        ExpMinusXVerifier.propagate_forward_exp_minus_x(getIntegrator());
    }

    @Test
    void propagate_forward_exp_to_cos_x() {
        propagate_exp_to_cos_x(Integrator.Direction.FORWARD);
    }

    @Test
    void propagate_backward_exp_to_cos_x() {
        propagate_exp_to_cos_x(Integrator.Direction.BACKWARD);
    }

    private void propagate_exp_to_cos_x(Integrator.Direction direction) {
        ExpCosXVerifier.propagate_exp_to_cos_x(getIntegrator(), direction, convergence_params_exp_to_cos_x());
    }

    // We assume errors are in the format |error| < (C/np)^a
    // The higher a , the better;
    // The smallest C, the better;
    abstract protected ConvergenceParams convergence_params_exp_to_cos_x();


}