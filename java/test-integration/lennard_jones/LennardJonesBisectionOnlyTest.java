package lennard_jones;

import schrodinger.solver.RefinementStrategy;

public class LennardJonesBisectionOnlyTest extends LennardJonesAbstractTest {

    public LennardJonesBisectionOnlyTest() {
        super(RefinementStrategy.BISECTION_ONLY, false);
        loadReferenceEnergies("resources/energies_lennard_jones_bisection.tsv");
    }

}
