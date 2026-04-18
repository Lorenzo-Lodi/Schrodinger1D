package lennard_jones;

import schrodinger.solver.RefinementStrategy;

public class LennardJonesBisectionThenRegulaFalsiTest extends LennardJonesAbstractTest {

    public LennardJonesBisectionThenRegulaFalsiTest() {
        super(RefinementStrategy.BISECTION_THEN_REGULA_FALSI, true);
    }

}
