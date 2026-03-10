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
        assertEquals(5.32756630584979, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(4.4120103459770625, val, 1e-14);
    }

}
