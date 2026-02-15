package schrodinger;

import org.junit.jupiter.api.Test;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.TaylorThreePoints;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialHarmonic;
import schrodinger.potential.SchrodingerSystem;
import schrodinger.solver.ShootingSolver;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class Test001 {

    @Test
    void test001() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(20, 1);
        Integrator integrator = new TaylorThreePoints();
        double mass = 2.0d;
        int nOfDesiredNodes = 0;
        int nOfPoints = 200;
        Grid grid = GridFactory.generateUniformGrid(14.0d, 26.0d, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, integrator);
        QuantumState ek = finder.findEigenvalue(nOfDesiredNodes);
        double refEnergy = 0.4997726295700777;
        assertEquals(ek.energy, refEnergy, 1e-12);
        assertTrue(ek.nodesUpper > refEnergy);
        assertTrue(ek.nodesLower < refEnergy);

        double refPerturbative = 2.270604341323613E-4;
        assertEquals(refPerturbative, ek.perturbativeCorrectionToEnergy, 1e-12);
    }

    @Test
    void test002() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(20, 1);
        Integrator integrator = new TaylorThreePoints();
        double mass = 2.0d;
        int nOfDesiredNodes = 10;
        int nOfPoints = 200;
        Grid grid = GridFactory.generateUniformGrid(14.0d, 26.0d, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, integrator);
        QuantumState ek = finder.findEigenvalue(nOfDesiredNodes);
        double refEnergy = 10.44953025724728;
        assertEquals(ek.energy, refEnergy, 1e-12);
        assertTrue(Math.abs(ek.energy - refEnergy) < 1e-12);
        assertTrue(ek.nodesUpper > refEnergy);
        assertTrue(ek.nodesLower < refEnergy);

        double refPerturbative = 0.050063577379135055;
        assertEquals(refPerturbative, ek.perturbativeCorrectionToEnergy, 1e-12);
    }

}
