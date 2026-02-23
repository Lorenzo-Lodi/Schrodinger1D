package schrodinger.integrator.expfitted;

import org.junit.jupiter.api.Test;
import schrodinger.integrator.AbstractIntegratorTest;
import schrodinger.integrator.Integrator;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test class for the ExponentiallyFitted integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class ExponentiallyFittedTest extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new ExponentiallyFitted();
    }

    @Test
    public void integrateForwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.FORWARD);
        assertEquals(1.80922655958763, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(1.8054002416708992, val, 1e-14);
    }

    @Test
    public void integrateForwardNegativeQTest() {
        double val = integrateOneStepNegativeQ(Integrator.Direction.FORWARD);
        assertEquals(-13.568451076776498, val, 1e-14);
    }

    @Test
    public void integrateBackwardNegativeQTest() {
        double val = integrateOneStepNegativeQ(Integrator.Direction.BACKWARD);
        assertEquals(-10.793967160597104, val, 1e-14);
    }

    @Test
    public void integrateForwardSmallZTest() {
        double val = integrateOneStepSmallZ(Integrator.Direction.FORWARD);
        assertEquals(-0.766890839498507, val, 1e-14);
    }

    @Test
    public void integrateBackwardSmallZTest() {
        double val = integrateOneStepSmallZ(Integrator.Direction.BACKWARD);
        assertEquals(-1.6824522709061138, val, 1e-14);
    }

}
