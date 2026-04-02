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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static schrodinger.PhysicalConstants.*;

public abstract class LennardJonesBisectionThenSecantTest {

    private static final double wellDepthInverseCm = 4050;
    private static final double wellDepthHartree = wellDepthInverseCm / HARTREE_TO_INVERSE_CM;
    private static final double rMinAng = 1.;
    private static final double rMinBohr = rMinAng / BOHR_TO_ANG;
    private static final PhysicalPotential potential = new PhysicalPotentialLennardJones(rMinBohr, wellDepthHartree, 6);
    private static final double mass = 16.85762920 * UMA_TO_ELECTRON_MASS;
    private static final RefinementStrategy strategy = RefinementStrategy.BISECTION_THEN_SECANT;
    private static final Map<Integer, Double> refEnergies = new HashMap<>();

    static {
        refEnergies.put(0, 371.40006854105536);
        refEnergies.put(1, 1046.35039758025960);
        refEnergies.put(2, 1636.86240202073500);
        refEnergies.put(3, 2147.54058595915200);
        refEnergies.put(4, 2583.19756725932500);
        refEnergies.put(5, 2948.86947887610630);
        refEnergies.put(6, 3249.83161605343780);
        refEnergies.put(7, 3491.61374341683600);
        refEnergies.put(8, 3680.01417388285970);
        refEnergies.put(9, 3821.11132786356030);
        refEnergies.put(10, 3921.27096821614170);
        refEnergies.put(11, 3987.14669051490550);
        refEnergies.put(12, 4025.67057896758570);
        refEnergies.put(13, 4044.03031326710800);
        refEnergies.put(14, 4049.62843497018100);
    }

    abstract Integrator getIntegrator();

    @Test
    void test_4000pts_00_nodes() {
        testCore(4000, 0);
    }

    @Test
    void test_4000pts_01_nodes() {
        testCore(4000, 1);
    }

    @Test
    void test_4000pts_02_nodes() {
        testCore(4000, 2);
    }

    @Test
    void test_4000pts_03_nodes() {
        testCore(4000, 3);
    }

    @Test
    void test_4000pts_04_nodes() {
        testCore(4000, 4);
    }

    @Test
    void test_4000pts_05_nodes() {
        testCore(4000, 5);
    }

    @Test
    void test_4000pts_06_nodes() {
        testCore(4000, 6);
    }

    @Test
    void test_4000pts_07_nodes() {
        testCore(4000, 7);
    }

    @Test
    void test_4000pts_08_nodes() {
        testCore(4000, 8);
    }

    @Test
    void test_4000pts_09_nodes() {
        testCore(4000, 9);
    }

    @Test
    void test_4000pts_10_nodes() {
        testCore(4000, 10);
    }

    @Test
    void test_4000pts_11_nodes() {
        testCore(4000, 11);
    }

    @Test
    void test_4000pts_12_nodes() {
        testCore(4000, 12);
    }

    @Test
    void test_4000pts_13_nodes() {
        testCore(4000, 13);
    }

    @Test
    void test_4000pts_14_nodes() {
        testCore(4000, 14);
    }

    private void testCore(int nOfPoints, int nOfDesiredNodes) {
        String className = this.getClass().getSimpleName().replaceAll("Test", "");
        OutputManager.initCommonOutputFile("LennardJones_" + className + "_" +
                nOfPoints + "_" + strategy.toString().toLowerCase() + ".log");
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
        assertEquals(refEnergies.get(nOfDesiredNodes), toInverseCm(ek.energy), 1e-3);

    }
}
