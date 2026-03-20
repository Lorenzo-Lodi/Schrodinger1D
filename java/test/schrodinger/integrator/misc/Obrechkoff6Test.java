package schrodinger.integrator.misc;

import org.junit.jupiter.api.Test;
import schrodinger.integrator.AbstractIntegratorTest;
import schrodinger.integrator.verification.ConvergenceParams;
import schrodinger.integrator.Integrator;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test class for the TaylorThreePoints integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class Obrechkoff6Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new Obrechkoff6();
    }


    @Test
    public void integrateForwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.FORWARD);
        assertEquals(2.3229953427486745, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(1.4410273443477517, val, 1e-14);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x() {
        return new ConvergenceParams(6.07, 0.70);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_minus_x() {
        return new ConvergenceParams(6.00, 2.1);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x_minus_x() {
        return new ConvergenceParams(6.00, 0.7);
    }

    @Override
    protected ConvergenceParams convergence_params_tanh_x() {
        return new ConvergenceParams(6.00, 1.8);
    }
}
