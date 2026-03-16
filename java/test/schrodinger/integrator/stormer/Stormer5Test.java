package schrodinger.integrator.stormer;

import org.junit.jupiter.api.Test;
import schrodinger.integrator.AbstractIntegratorTest;
import schrodinger.integrator.verification.ConvergenceParams;
import schrodinger.integrator.Integrator;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test class for the Stormer5 integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class Stormer5Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new Stormer5();
    }

    @Test
    public void integrateForwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.FORWARD);
        assertEquals(2.6789537448608023, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(2.3756666670094093, val, 1e-14);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x() {
        return new ConvergenceParams(4.85, 28.6);
    }

}
