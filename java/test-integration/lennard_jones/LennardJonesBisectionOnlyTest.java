package lennard_jones;

import schrodinger.integrator.IntegratorFactory;
import schrodinger.solver.RefinementStrategy;

public class LennardJonesBisectionOnlyTest extends LennardJonesAbstractTest {

    public LennardJonesBisectionOnlyTest() {
        super(RefinementStrategy.BISECTION_ONLY, false);
    }

}
