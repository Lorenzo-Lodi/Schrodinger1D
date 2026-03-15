package schrodinger.integrator;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class RK45DPTest extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new RK45DP();
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

    @Override
    protected double error_bound_propagate_exp_to_cos_x(int np) {
        double  a = 5.95;
        double C = 14.3;
        return Math.pow(C / np, a);
    }

}
