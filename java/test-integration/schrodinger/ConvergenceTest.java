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

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class ConvergenceTest {

    @Test
    void eigenvalue_errors_for_harmonic_energy() {
        double mass = 2.0d;
        PhysicalPotential potential = new PhysicalPotentialHarmonic(20, 1);
        int nOfDesidedNodes = 15;
        List<Integrator> integrators = IntegratorFactory.getAll();
        integrators = integrators.stream().filter((x) -> x.minHistoryLength() <= 2).collect(Collectors.toList());

        for (Integrator integrator : integrators) {
            System.out.print(integrator.getClass().getSimpleName() + " ");
        }
        System.out.println();

        for (int nOfPoints = 40; nOfPoints <= 500; nOfPoints += 10) {
            Grid grid = GridFactory.generateUniformGrid(14.0d, 26.0d, nOfPoints);
            SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
            List<ShootingSolver> solvers = new ArrayList<>();
            for (Integrator integrator : integrators) {
                solvers.add(new ShootingSolver(system, integrator));
            }
            double exactEnergy = 0.5 + nOfDesidedNodes;
            List<Double> energies = new ArrayList<>();
            for (ShootingSolver s : solvers) {
                double energy = s.findEigenvalue(nOfDesidedNodes).energy;
                energies.add(energy);
            }

            System.out.print(padInt(nOfPoints) + " ");
            for (Double error : energies) {
                System.out.print(padFloat(error) + " ");
            }
            System.out.println();
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
