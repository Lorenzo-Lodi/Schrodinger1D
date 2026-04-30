package lennard_jones;

import schrodinger.solver.RefinementStrategy;

public class LennardJonesBisectionThenAndersonBjorckTest extends LennardJonesAbstractTest {

    public LennardJonesBisectionThenAndersonBjorckTest() {
        super(RefinementStrategy.BISECTION_THEN_ANDERSON_BJORCK, false);
    }

}
