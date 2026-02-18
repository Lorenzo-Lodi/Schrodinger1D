package schrodinger;

import org.junit.jupiter.api.Test;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.IntegratorFactory;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialHarmonic;
import schrodinger.potential.SchrodingerSystem;
import schrodinger.solver.ShootingSolver;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class ConvergenceTest {

    @Test
    void test001() {
        double mass = 2.0d;
//        int nOfPoints = 200;
        PhysicalPotential potential = new PhysicalPotentialHarmonic(20, 1);
        int nOfDesidedNodes = 15;

        for (int nOfPoints = 100; nOfPoints <= 1000; nOfPoints += +100) {
            Grid grid = GridFactory.generateUniformGrid(14.0d, 26.0d, nOfPoints);
            SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
            ShootingSolver finder1 = new ShootingSolver(system, IntegratorFactory.getTaylorThreePoints());
            ShootingSolver finder2 = new ShootingSolver(system, IntegratorFactory.getNumerov());
            ShootingSolver finder3 = new ShootingSolver(system, IntegratorFactory.getVignoli());
            ShootingSolver finder4 = new ShootingSolver(system, IntegratorFactory.getExponentiallyFitted());

            QuantumState e1 = finder1.findEigenvalue(nOfDesidedNodes);
            QuantumState e2 = finder2.findEigenvalue(nOfDesidedNodes);
            QuantumState e3 = finder3.findEigenvalue(nOfDesidedNodes);
            QuantumState e4 = finder4.findEigenvalue(nOfDesidedNodes);

            double exact = 0.5 + nOfDesidedNodes;
            double err1 = exact - e1.energy;
            double err1p = exact - e1.energy - e1.perturbativeCorrectionToEnergy;
            double err2 = exact - e2.energy;
            double err2p = exact - e2.energy - e2.perturbativeCorrectionToEnergy;
            double err3 = exact - e3.energy;
            double err4 = exact - e4.energy;
            System.out.println(padInt(nOfPoints) + " " + padFloat(err1) + " " + padFloat(err1p)
                    + " " + padFloat(err2) + " " + padFloat(err2p)
                    + " " + padFloat(err3)
                    + " " + padFloat(err4)
            );
        }

    }

    private static String padFloat(double number) {
        return padFloat(number, 15, 18);
    }

    private static String padInt(int number) {
        return padInt(number, 6);
    }

    private static String padding(double number, int width) {
        int magnitude = number <= 1 ? 0 : (int) Math.log10(Math.abs(number));
        String padding = "";
        int signPadding = (number >= 0) ? 0 : 1;
        for (int i = 0; i < width - magnitude - signPadding; i++) {
            padding += " ";
        }
        return padding;
    }

    private static String padFloat(Double number, int nOfDecimals, int width) {
        return padding(number, width - nOfDecimals - 2) + String.format("%." + nOfDecimals + "f", number);
    }

    private static String padInt(Integer number, int width) {
        return padding(number, width) + number;
    }

}
