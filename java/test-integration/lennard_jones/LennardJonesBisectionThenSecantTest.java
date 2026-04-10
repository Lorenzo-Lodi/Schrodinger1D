package lennard_jones;

import schrodinger.integrator.Integrator;
import schrodinger.integrator.IntegratorFactory;
import schrodinger.solver.RefinementStrategy;

import java.util.List;

public class LennardJonesBisectionThenSecantTest extends LennardJonesAbstractTest {

    public LennardJonesBisectionThenSecantTest() {
        super(RefinementStrategy.BISECTION_THEN_SECANT, false);
    }

}
