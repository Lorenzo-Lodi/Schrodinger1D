package lennard_jones;

import schrodinger.integrator.Integrator;
import schrodinger.integrator.numerovlike.EFNFixedBeta;

public class EFNFixedBetaTest extends LennardJonesBisectionOnlyTest {

    @Override
    Integrator getIntegrator() {
        return new EFNFixedBeta();
    }
}
