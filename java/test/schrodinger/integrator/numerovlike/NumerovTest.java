package schrodinger.integrator.numerovlike;

import org.junit.jupiter.api.Test;
import schrodinger.integrator.AbstractIntegratorTest;
import schrodinger.integrator.ConvergenceParams;
import schrodinger.integrator.Integrator;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test class for the Numerov integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class NumerovTest extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new Numerov();
    }

    @Test
    public void integrateForwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.FORWARD);
        assertEquals(2.620986993887933, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(2.3476065744871613, val, 1e-14);
    }

    @Test
    public void integrateForwardNegativeQTest() {
        double val = integrateOneStepNegativeQ(Integrator.Direction.FORWARD);
        assertEquals(-31.235813595724405, val, 1e-14);
    }

    @Test
    public void integrateBackwardNegativeQTest() {
        double val = integrateOneStepNegativeQ(Integrator.Direction.BACKWARD);
        assertEquals(-12.180183007186118, val, 1e-14);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x() {
        return new ConvergenceParams(3.62, 15.1);
    }


}
