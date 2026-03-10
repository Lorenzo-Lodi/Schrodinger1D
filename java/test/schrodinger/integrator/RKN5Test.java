package schrodinger.integrator;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class RKN5Test extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new RKN5();
    }

//    @Test
//    public void integrateForwardReferenceTest() {
//        double val = integrateOneStep(Integrator.Direction.FORWARD);
//        assertEquals(0.36588338422771727, val, 1e-14);
//    }
//
//    @Test
//    public void integrateBackwardReferenceTest() {
//        double val = integrateOneStep(Integrator.Direction.BACKWARD);
//        assertEquals(0.34054956612279075, val, 1e-14);
//    }

}
