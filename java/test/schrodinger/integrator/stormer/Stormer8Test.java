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
public class Stormer8Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new Stormer8();
    }
    @Test
    public void integrateForwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.FORWARD);
        assertEquals(4.718050189782055, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(4.444088283345244, val, 1e-14);
    }

    @Override
    protected double error_bound_propagate_exp_to_cos_x(int np) {
        double a = 8.63;
        double C = 52.1;
        return Math.pow(C / np, a);
    }

}
