package schrodinger;

import org.junit.jupiter.api.Test;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.IntegratorFactory;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialHarmonic;
import schrodinger.potential.SchrodingerSystem;
import schrodinger.integrator.pt_correction.NumerovPTCorrector;
import schrodinger.integrator.pt_correction.PTCorrector;
import schrodinger.integrator.pt_correction.TaylorThreePointsPTCorrector;
import schrodinger.solver.ShootingSolver;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class IntegrationTest {

    @Test
    void test001() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(20, 1);
        Integrator integrator = IntegratorFactory.getTaylorThreePoints();
        double mass = 2.0d;
        int nOfDesiredNodes = 0;
        int nOfPoints = 200;
        Grid grid = GridFactory.generateUniformGrid(14.0d, 26.0d, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, integrator);
        QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes);
        double refEnergy = 0.4997726295700777;
        assertEquals(refEnergy, ek.energy, 1e-12);
        assertTrue(ek.nodesUpper > refEnergy);
        assertTrue(ek.nodesLower < refEnergy);

        double refPerturbative = 2.270604341323613E-4;
        PTCorrector pt = new TaylorThreePointsPTCorrector();
        pt.computeAndSet(ek);
        assertEquals(refPerturbative, ek.perturbativeCorrectionToEnergy, 1e-12);
    }

    @Test
    void test002() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(20, 1);
        Integrator integrator = IntegratorFactory.getTaylorThreePoints();
        double mass = 2.0d;
        int nOfDesiredNodes = 10;
        int nOfPoints = 200;
        Grid grid = GridFactory.generateUniformGrid(14.0d, 26.0d, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, integrator);
        QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes);
        double refEnergy = 10.44953025724728;
        assertEquals(refEnergy, ek.energy, 1e-12);
        assertTrue(Math.abs(ek.energy - refEnergy) < 1e-12);
        assertTrue(ek.nodesUpper > refEnergy);
        assertTrue(ek.nodesLower < refEnergy);

        double refPerturbative = 0.050063577379135055;
        PTCorrector pt = new TaylorThreePointsPTCorrector();
        pt.computeAndSet(ek);
        assertEquals(refPerturbative, ek.perturbativeCorrectionToEnergy, 1e-12);
    }


    @Test
    void test003() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(20, 1);
        Integrator integrator = IntegratorFactory.getNumerov();
        double mass = 2.0d;
        int nOfDesiredNodes = 0;
        int nOfPoints = 200;
        Grid grid = GridFactory.generateUniformGrid(14.0d, 26.0d, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, integrator);
        QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes);
        double refEnergy = 0.49999979318999277;
        assertEquals(refEnergy, ek.energy, 1e-12);
        assertTrue(ek.nodesUpper > refEnergy);
        assertTrue(ek.nodesLower < refEnergy);

        double refPerturbative = 2.0485571986324924E-7;
        PTCorrector pt = new NumerovPTCorrector();
        pt.computeAndSet(ek);
        assertEquals(refPerturbative, ek.perturbativeCorrectionToEnergy, 1e-12);
    }

    @Test
    void test004() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(20, 1);
        Integrator integrator = IntegratorFactory.getNumerov();
        double mass = 2.0d;
        int nOfDesiredNodes = 10;
        int nOfPoints = 200;
        Grid grid = GridFactory.generateUniformGrid(14.0d, 26.0d, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, integrator);
        QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes);
        double refEnergy = 10.499675728015717;
        assertEquals(refEnergy, ek.energy, 1e-12);
        assertTrue(Math.abs(ek.energy - refEnergy) < 1e-12);
        assertTrue(ek.nodesUpper > refEnergy);
        assertTrue(ek.nodesLower < refEnergy);

        double refPerturbative = 3.0814906265290394E-4;
        PTCorrector pt = new NumerovPTCorrector();
        pt.computeAndSet(ek);
        assertEquals(refPerturbative, ek.perturbativeCorrectionToEnergy, 1e-12);
    }

    @Test
    void test005() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(20, 1);
        Integrator integrator = IntegratorFactory.getEFNFixedBeta();
        double mass = 2.0d;
        int nOfDesiredNodes = 0;
        int nOfPoints = 200;
        Grid grid = GridFactory.generateUniformGrid(14.0d, 26.0d, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, integrator);
        QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes);
        double refEnergy = 0.4999997794980161;
        assertEquals(refEnergy, ek.energy, 1e-12);
        assertTrue(ek.nodesUpper > refEnergy);
        assertTrue(ek.nodesLower < refEnergy);

        double refPerturbative = 0.; // Unavailable for this method.
        assertEquals(refPerturbative, ek.perturbativeCorrectionToEnergy, 1e-12);
    }

    @Test
    void test006() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(20, 1);
        Integrator integrator = IntegratorFactory.getEFNFixedBeta();
        double mass = 2.0d;
        int nOfDesiredNodes = 10;
        int nOfPoints = 200;
        Grid grid = GridFactory.generateUniformGrid(14.0d, 26.0d, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, integrator);
        QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes);
        double refEnergy = 10.499995315122833;
        assertEquals(refEnergy, ek.energy, 1e-12);
        assertTrue(Math.abs(ek.energy - refEnergy) < 1e-12);
        assertTrue(ek.nodesUpper > refEnergy);
        assertTrue(ek.nodesLower < refEnergy);

        double refPerturbative = 0;  // Unavailable for this method.
        assertEquals(refPerturbative, ek.perturbativeCorrectionToEnergy, 1e-12);
    }

    @Test
    void test007() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(20, 1);
        Integrator integrator = IntegratorFactory.getNumerov();
        double mass = 2.0d;
        int nOfDesiredNodes = 10;
        int nOfPoints = 200;
        Grid grid = GridFactory.generateLogarithmicGrid(14.0d, 26.0d, nOfPoints, 20.);
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, integrator);
        QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes);
        double refEnergy = 10.499663998737587;
        assertEquals(refEnergy, ek.energy, 1e-12);
        assertTrue(Math.abs(ek.energy - refEnergy) < 1e-12);
        assertTrue(ek.nodesUpper > refEnergy);
        assertTrue(ek.nodesLower < refEnergy);

        double refPerturbative = 3.1895686240242146E-4;
        PTCorrector pt = new NumerovPTCorrector();
        pt.computeAndSet(ek);
        assertEquals(refPerturbative, ek.perturbativeCorrectionToEnergy, 1e-12);
    }

    @Test
    void test008() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(20, 1);
        Integrator integrator = IntegratorFactory.getNumerov();
        double mass = 2.0d;
        int nOfDesiredNodes = 10;
        int nOfPoints = 200;
        Grid grid = GridFactory.generateSqrtGrid(14.0d, 26.0d, nOfPoints, 20.);
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, integrator);
        QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes);
        double refEnergy = 10.499649544521517;
        assertEquals(refEnergy, ek.energy, 1e-12);
        assertTrue(Math.abs(ek.energy - refEnergy) < 1e-12);
        assertTrue(ek.nodesUpper > refEnergy);
        assertTrue(ek.nodesLower < refEnergy);

        double refPerturbative = 3.11832250753927E-4;
        PTCorrector pt = new NumerovPTCorrector();
        pt.computeAndSet(ek);
        assertEquals(refPerturbative, ek.perturbativeCorrectionToEnergy, 1e-12);
    }

    @Test
    void test009() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(20, 1);
        Integrator integrator = IntegratorFactory.getNumerov();
        double mass = 2.0d;
        int nOfDesiredNodes = 10;
        int nOfPoints = 200;
        Grid grid = GridFactory.generateSurkusGrid(14.0d, 26.0d, nOfPoints, 20., 1.1);
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        ShootingSolver finder = new ShootingSolver(system, integrator);
        QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes);
        double refEnergy = 10.499638433533447;
        assertEquals(refEnergy, ek.energy, 1e-12);
        assertTrue(Math.abs(ek.energy - refEnergy) < 1e-12);
        assertTrue(ek.nodesUpper > refEnergy);
        assertTrue(ek.nodesLower < refEnergy);

        double refPerturbative = 3.424279849395494E-4;
        PTCorrector pt = new NumerovPTCorrector();
        pt.computeAndSet(ek);
        assertEquals(refPerturbative, ek.perturbativeCorrectionToEnergy, 1e-12);
    }

}
