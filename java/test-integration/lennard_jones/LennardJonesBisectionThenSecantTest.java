package lennard_jones;

import schrodinger.solver.RefinementStrategy;

public class LennardJonesBisectionThenSecantTest extends LennardJonesAbstractTest {

    public LennardJonesBisectionThenSecantTest() {
        super(RefinementStrategy.BISECTION_THEN_SECANT, true);
    }

}
