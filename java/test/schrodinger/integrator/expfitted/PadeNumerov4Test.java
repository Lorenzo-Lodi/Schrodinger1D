package schrodinger.integrator.expfitted;

import schrodinger.integrator.AbstractIntegratorTest;
import schrodinger.integrator.Integrator;

public class PadeNumerov4Test extends AbstractIntegratorTest {
    @Override
    protected Integrator getIntegrator() {
        return new PadeNumerov4();
    }
}
