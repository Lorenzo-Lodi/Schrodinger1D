package schrodinger.integrator;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test class for the TaylorThreePoints integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class TaylorThreePointsTest extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new TaylorThreePoints();
    }

    @Test
    public void integrateForwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.FORWARD);
        assertEquals(5.32756630584979, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(4.4120103459770625, val, 1e-14);
    }

    @Override
    protected double error_bound_propagate_exp_to_cos_x(int np) {
        double a = 1.88;
        double C = 85.4;
        return Math.pow(C / np, a);
    }

}
