package schrodinger.solver;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import schrodinger.SchrodingerSystem;
import schrodinger.grid.Grid;
import schrodinger.grid.MappingStrategy;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.IntegratorFactory;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialMorse;


public class ShootingSolverTest {
    Grid grid;
    PhysicalPotential potential;
    SchrodingerSystem system;
    Integrator integrator;
    ShootingSolver solver;

    @BeforeEach
    public void setup() {
        grid = new Grid(0.5, 2.5, 20, new MappingStrategy() {
        });
        potential = new PhysicalPotentialMorse(1.1, 0.8, 0.1);
        system = new SchrodingerSystem(potential, 2, grid);
        integrator = IntegratorFactory.getNumerov();
        solver = new ShootingSolver(system, integrator);
    }

    @Test
    public void test01() {

    }

}
