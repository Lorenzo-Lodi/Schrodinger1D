package schrodinger.integrator;

public class Cowell4Test extends AbstractIntegratorTest {
    @Override
    protected Integrator getIntegrator() {
        return new Cowell4();
    }
}
