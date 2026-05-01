package morse;

import schrodinger.OutputManager;
import schrodinger.QuantumLevel;
import schrodinger.SchrodingerSystem;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.IntegratorFactory;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialLennardJones;
import schrodinger.potential.PhysicalPotentialMorse;
import schrodinger.solver.RefinementStrategy;
import schrodinger.solver.ShootingSolver;

import java.util.HashMap;
import java.util.Map;

import static schrodinger.PhysicalConstants.*;

public class MorseMain {

    public static void main(String[] args) {

        PhysicalPotential potential = new PhysicalPotentialMorse(2.0, 1.5, toHartree(50000.));
        double xmin = 1.0;
        double xmax = 12.;
        int nOfPoints = (int) ((xmax - xmin) / 0.005);
        System.out.println(toInverseCm(potential.value(xmin))  + " " + toInverseCm(potential.value(xmax))  + "  " + nOfPoints);

        Grid grid = GridFactory.generateUniformGrid(xmin, xmax, nOfPoints);
        double mass = 100. * UMA_TO_ELECTRON_MASS;
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);

        Integrator integrator = IntegratorFactory.getNumerov();
        ShootingSolver finder = new ShootingSolver(system, integrator);

        for (int nOfDesiredNodes = 0; nOfDesiredNodes <= 200; nOfDesiredNodes++) {
            QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes, RefinementStrategy.BISECTION_ONLY);
            System.out.println(nOfDesiredNodes + " " + toInverseCm(ek.energy));
        }


    }

}
