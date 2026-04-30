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
public class PredictorCorrector8i1Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new PredictorCorrector8i1();
    }

    @Test
    public void integrateForwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.FORWARD);
        assertEquals(1.7720924944425227, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(1.9871191175009, val, 1e-14);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x() {
        return new ConvergenceParams(8.00, 0.45);
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
        return new ConvergenceParams(8.00, 0.8);
    }

}
