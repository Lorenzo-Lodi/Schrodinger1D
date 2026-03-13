package schrodinger.integrator;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class RKN4Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new RKN4();
    }

    @Test
    public void integrateForwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.FORWARD);
        assertEquals(0.3996918280374979, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(0.336941647798199, val, 1e-14);
    }

}
