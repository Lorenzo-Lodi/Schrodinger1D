package schrodinger;

import org.junit.jupiter.api.Test;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.notworking.Cowell8;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.IntegratorFactory;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialHarmonic;
import schrodinger.potential.SchrodingerSystem;
import schrodinger.solver.RefinementStrategy;
import schrodinger.solver.ShootingSolver;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class NewIntegratorTest {

    @Test
    void test001() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(20, 1);
        Integrator integrator = IntegratorFactory.getExponentiallyFitted();
        double mass = 2.0d;
        int nOfDesiredNodes = 0;
        int nOfPoints = 200;
        Grid grid = GridFactory.generateUniformGrid(14.0d, 26.0d, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, integrator);
        QuantumState ek = finder.findEigenvalue(nOfDesiredNodes);
        double refEnergy = 0.4999997794980129;
        assertEquals(refEnergy, ek.energy, 1e-12);
        assertTrue(ek.nodesUpper > refEnergy);
        assertTrue(ek.nodesLower < refEnergy);

        double refPerturbative = 0.;
        assertEquals(refPerturbative, ek.perturbativeCorrectionToEnergy, 1e-12);
    }

    @Test
    void test002() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(20, 1);
        Integrator integrator = new Cowell8();
        double mass = 2.0d;
        int nOfDesiredNodes = 0;
        int nOfPoints = 200;
        Grid grid = GridFactory.generateUniformGrid(14.0d, 26.0d, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, integrator);
        QuantumState ek = finder.findEigenvalue(nOfDesiredNodes, RefinementStrategy.BISECTION_ONLY);
        System.out.println(ek.energy);
//        double refEnergy = 0.4999997794980129;
//        assertEquals(refEnergy, ek.energy, 1e-12);
//        assertTrue(ek.nodesUpper > refEnergy);
//        assertTrue(ek.nodesLower < refEnergy);
//
//        double refPerturbative = 0.;
//        assertEquals(refPerturbative, ek.perturbativeCorrectionToEnergy, 1e-12);
    }


}
