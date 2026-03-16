package schrodinger.integrator;

import org.junit.jupiter.api.Test;
import schrodinger.integrator.rungekutta.RK45DP;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class RK45DPTest extends AbstractIntegratorTest {

    @Override
    protected Integrator getIntegrator() {
        return new RK45DP();
    }

    @Test
    public void integrateForwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.FORWARD);
        assertEquals(0.9418954134438314, val, 1e-14);
    }

    @Test
    public void integrateBackwardReferenceTest() {
        double val = integrateOneStep(Integrator.Direction.BACKWARD);
        assertEquals(0.8376003369215926, val, 1e-14);
    }

    @Override
    protected ConvergenceParams convergence_params_exp_to_cos_x() {
        return new ConvergenceParams(5.95, 14.3);
    }

}
