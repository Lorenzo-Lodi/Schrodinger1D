package schrodinger.integrator.predcorr;

import org.junit.jupiter.api.Test;
import schrodinger.integrator.AbstractIntegratorTest;
import schrodinger.integrator.Integrator;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test class for the PredictorCorrector6 integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class PredictorCorrector6Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new PredictorCorrector6();
    }

    @Test
    public void integrateForwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.FORWARD);
        assertEquals(2.183237656208287, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(2.105467768082526, val, 1e-14);
    }

    @Override
    protected double error_bound_propagate_exp_to_cos_x(int np) {
        double a = 8.88;
        double C = 41.2;
        return Math.pow(C / np, a);
    }

}