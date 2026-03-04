package schrodinger.integrator.predcorr;

import org.junit.jupiter.api.Test;
import schrodinger.integrator.AbstractIntegratorTest;
import schrodinger.integrator.Integrator;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Test class for the PredictorCorrector8Abstract integrator.
 * Tests the convergence of the integrator by comparing numerical solutions
 * with the exact analytical solution for a harmonic oscillator potential.
 */
public class PredictorCorrector8NumerovIter2Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new PredictorCorrector8NumerovIter2();
    }
    @Test
    public void integrateForwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.FORWARD);
        assertEquals(-1.8563504317777002, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(-2.6719729048822267, val, 1e-14);
    }
}
