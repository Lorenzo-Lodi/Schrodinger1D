package lennard_jones;

import schrodinger.solver.RefinementStrategy;

public class LennardJonesBisectionThenBidirectionalTest extends LennardJonesAbstractTest {

    public LennardJonesBisectionThenBidirectionalTest() {
        super(RefinementStrategy.BISECTION_THEN_BIDIRECTIONAL, false);
        loadReferenceEnergies("resources/energies_lennard_jones_improved_regula_falsi.tsv");
    }

}
