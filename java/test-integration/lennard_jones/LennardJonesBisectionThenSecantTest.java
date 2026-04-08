package lennard_jones;

import schrodinger.integrator.Integrator;
import schrodinger.integrator.IntegratorFactory;
import schrodinger.solver.RefinementStrategy;

import java.util.List;

public class LennardJonesBisectionThenSecantTest extends LennardJonesAbstractTest {

    public LennardJonesBisectionThenSecantTest() {
        super(RefinementStrategy.BISECTION_THEN_SECANT, true);

        // For now focus on Obrechkoff6, one of the more problematic, especially for large xmax
        this.integrators = List.of(IntegratorFactory.getObrechkoff6());
    }

}
