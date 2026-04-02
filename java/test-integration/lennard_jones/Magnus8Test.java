package lennard_jones;

import schrodinger.integrator.Integrator;
import schrodinger.integrator.IntegratorFactory;

public class Magnus8Test extends LennardJonesBisectionThenSecantTest {

    @Override
    Integrator getIntegrator() {
        return IntegratorFactory.getCFMagnus8();
    }
}
