package schrodinger.integrator.stormer;

import org.junit.jupiter.api.Test;
import schrodinger.integrator.AbstractIntegratorTest;
import schrodinger.integrator.Integrator;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test class for the Stormer6 integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class Stormer8iTest extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new Stormer8i();
    }
    @Test
    public void integrateForwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.FORWARD);
        assertEquals(2.819478335032403, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(2.527657332202554, val, 1e-14);
    }

    @Override
    protected double error_bound_propagate_exp_to_cos_x(int np) {
        double a = 8.70;
        double C = 36.6;
        return Math.pow(C / np, a);
    }

}
