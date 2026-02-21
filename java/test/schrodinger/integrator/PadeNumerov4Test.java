package schrodinger.integrator;

public class PadeNumerov4Test extends AbstractIntegratorTest {
    @Override
    protected Integrator getIntegrator() {
        return new PadeNumerov4();
    }
}
