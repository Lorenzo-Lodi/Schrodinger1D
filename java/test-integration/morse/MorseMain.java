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

        double a = 1.5;
        double De = toHartree(50000.);

        PhysicalPotential potential = new PhysicalPotentialMorse(2.0, a, De);
        double xmin = 1.0;
        double xmax = 12.;
        int nOfPoints = (int) ((xmax - xmin) / 0.005);

        Grid grid = GridFactory.generateUniformGrid(xmin, xmax, nOfPoints);
        double mass = 100. * UMA_TO_ELECTRON_MASS;
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);

        Integrator integrator = IntegratorFactory.getPC8i2();
        ShootingSolver finder = new ShootingSolver(system, integrator);

        double omega0 = a * Math.sqrt(2.0 * De / mass);
        double xe = omega0 / (4. * De);
        double A = 1. / xe;
        int nOfBoundStates = (int) ((A + 1.) * 0.5);
        System.out.printf("omega0=%20.8f cm-1; xe=%20.8f nOfBoundStates=%10d\n", toInverseCm(omega0), xe, nOfBoundStates);

        System.out.println(toInverseCm(potential.value(xmin)) + " " + toInverseCm(potential.value(xmax)) + "  " + nOfPoints);

        System.out.printf("%10s %20s %20s %20s \n", "n", "calc", "exact", "exact-calc");
        for (int nOfDesiredNodes = 0; nOfDesiredNodes <= 200; nOfDesiredNodes++) {
            QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes, RefinementStrategy.BISECTION_ONLY);
            double exact = omega0 * (nOfDesiredNodes + 0.5) * (1. - xe * (nOfDesiredNodes + 0.5));
            double diff = exact - ek.energy;
            System.out.printf("%10d %20.8f %20.8f %20.8f\n", nOfDesiredNodes, toInverseCm(ek.energy), toInverseCm(exact), toInverseCm(diff));
        }


    }

}
