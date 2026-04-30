package schrodinger.integrator.multi_step.predcorr;

import org.junit.jupiter.api.Test;
import schrodinger.integrator.AbstractIntegratorTest;
import schrodinger.integrator.verification.ConvergenceParams;
import schrodinger.integrator.Integrator;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test class for the PredictorCorrector8Abstract integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class PredictorCorrector8i2Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new PredictorCorrector8i2();
    }

    @Test
    public void integrateForwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.FORWARD);
        assertEquals(1.9597799044581992, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(2.0023768806193236, val, 1e-14);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x() {
        return new ConvergenceParams(6.33, 0.8);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_minus_x() {
        return new ConvergenceParams(8.00, 0.7);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x_minus_x() {
        return new ConvergenceParams(8.00, 0.5);
    }

    @Override
    protected ConvergenceParams convergence_params_tanh_x() {
        return new ConvergenceParams(8.00, 1.1);
    }
}
