package lennard_jones;

import org.junit.jupiter.api.Test;
import schrodinger.OutputManager;
import schrodinger.QuantumLevel;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.Integrator;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialLennardJones;
import schrodinger.potential.SchrodingerSystem;
import schrodinger.solver.RefinementStrategy;
import schrodinger.solver.ShootingSolver;

import java.util.HashMap;
import java.util.Map;

import static schrodinger.PhysicalConstants.*;

public abstract class LennardJonesBisectionOnlyTest {

    private static final double wellDepthInverseCm = 4050;
    private static final double wellDepthHartree = wellDepthInverseCm / HARTREE_TO_INVERSE_CM;
    private static final double rMinAng = 1.;
    private static final double rMinBohr = rMinAng / BOHR_TO_ANG;
    private static final PhysicalPotential potential = new PhysicalPotentialLennardJones(rMinBohr, wellDepthHartree, 6);
    private static final double mass = 16.85762920 * UMA_TO_ELECTRON_MASS;
    private static final RefinementStrategy strategy = RefinementStrategy.BISECTION_ONLY;
    private static final Map<Integer, Double> refEnergies = new HashMap<>();

    static {
        refEnergies.put(0, -3678.599931);
        refEnergies.put(1, -3003.649602);
        refEnergies.put(2, -2413.137598);
        refEnergies.put(3, -1902.459414);
        refEnergies.put(4, -1466.802433);
        refEnergies.put(5, -1101.130521);
        refEnergies.put(6, -800.168384);
        refEnergies.put(7, -558.386257);
        refEnergies.put(8, -369.985826);
        refEnergies.put(9, -228.888672);
        refEnergies.put(10, -128.729032);
        refEnergies.put(11, -62.853309);
        refEnergies.put(12, -24.329421);
        refEnergies.put(13, -5.969687);
        refEnergies.put(14, -0.371565);
    }

    abstract Integrator getIntegrator();

    @Test
    void test_300pts_0_nodes() {
        int nOfPoints = 300;
        int nOfDesiredNodes = 0;
        OutputManager.initCommonOutputFile("DissociationTest_" + this.getClass().getSimpleName()
                + nOfPoints + "_" + strategy.toString().toLowerCase() + ".log");
        Grid grid = GridFactory.generateUniformGrid(1.2d, 45., nOfPoints);
        OutputManager.write("Grid info");
        OutputManager.write(String.format("%10s %22s %22s %22s", "i", "r", "y", "V"));
        for (int i = 0; i < grid.getNumberOfPoints(); i++) {
            double v = potential.value(grid.rAtGridPoint(i));
            OutputManager.write(String.format("%10d %22.8f %22.8f %25.6f",
                    i, grid.rAtGridPoint(i), grid.yAtGridPoint(i), toInverseCm(v)));
        }
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, getIntegrator());

        QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes, strategy);

    }
}
