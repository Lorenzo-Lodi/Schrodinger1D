package lennard_jones;

import schrodinger.solver.RefinementStrategy;

public class LennardJonesBisectionThenAndersonBjorckTest extends LennardJonesAbstractTest {

    public LennardJonesBisectionThenAndersonBjorckTest() {
        super(RefinementStrategy.BISECTION_THEN_ANDERSON_BJORCK, false);
        loadReferenceEnergies("resources/energies_lennard_jones_improved_regula_falsi.tsv");
    }

}
