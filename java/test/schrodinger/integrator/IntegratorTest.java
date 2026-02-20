package schrodinger.integrator;

import org.junit.jupiter.api.Test;
import schrodinger.QuantumState;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialHarmonic;
import schrodinger.potential.SchrodingerSystem;

import java.util.ArrayList;
import java.util.List;

public class IntegratorTest {

    @Test
    public void testIntegrator() {
//        Integrator integrator = new TaylorThreePoints(); // CHECKED: order 2
//        Integrator integrator = new Numerov(); // CHECKED: order 4
//        Integrator integrator = new Vignoli(); // CHECKED: order 4
//        Integrator integrator = new RaptisAllison(); // CHECKED: order 4
//        Integrator integrator = new ExponentiallyFitted(); // CHECKED: order 4

//        Integrator integrator = new PredictorCorrector6(); // NEEDS ONE FURTHER POINTS! Order about 6

//        Integrator integrator = new PredictorCorrector8(); // NEEDS TWO FURTHER POINTS! Order about 6
        Integrator integrator = new PredictorCorrector8Alt(); // NEEDS TWO FURTHER POINTS!


//        Integrator integrator = new Stormer7(); // NEEDS TWO FURTHER POINTS! Order about 4.8 or so!

//        Integrator integrator = new Stormer8(); // NEED THREE FURTHER POINTS! Order about 5.8 or so!

        tracking_harmonic_ground_state(integrator);
    }

    public static void tracking_harmonic_ground_state(Integrator integrator) {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(10, 1.0);

        List<Double> errors_points_over_10 = new ArrayList<>();
        List<Double> errors_points_over_4 = new ArrayList<>();
        List<Double> errors_points_over_2 = new ArrayList<>();

        for (int nOfPoints = 100; nOfPoints <= 1000; nOfPoints += 100) {
            Grid grid = GridFactory.generateUniformGrid(6.0, 14.0, nOfPoints);
            SchrodingerSystem system = new SchrodingerSystem(potential, 2., grid);
            QuantumState state = new QuantumState(system);
            state.energy = 0.5;
            state.psi[0] = exactSolution(grid.getRValue(0));
            state.psi[1] = exactSolution(grid.getRValue(1));

            state.psi[2] = exactSolution(grid.getRValue(2));
            state.psi[3] = exactSolution(grid.getRValue(3));
//            state.psi[4] = exactSolution(grid.getRValue(4));

            for (int n = 3; n < nOfPoints - 1; n++) {
                state.psi[n + 1] = integrator.propagate(state.psi, n, state, Integrator.Direction.FORWARD);
            }

            errors_points_over_10.add(exactSolution(grid.getRValue(nOfPoints / 10)) - state.psi[nOfPoints / 10]);
            errors_points_over_4.add(exactSolution(grid.getRValue(nOfPoints / 4)) - state.psi[nOfPoints / 4]);
            errors_points_over_2.add(exactSolution(grid.getRValue(nOfPoints / 2)) - state.psi[nOfPoints / 2]);

        }

        int i = 0;
        for (int nOfPoints = 100; nOfPoints <= 1000; nOfPoints += 100) {
            System.out.println(nOfPoints + " " + errors_points_over_10.get(i) + " " + errors_points_over_4.get(i) + " " + errors_points_over_2.get(i));
            i++;
        }
        System.out.println();
//        int i = 1;
        i=1;
        for (int nOfPoints = 200; nOfPoints <= 1000; nOfPoints += 100) {
            double denom = Math.log10(((double) nOfPoints) / ((double) nOfPoints - 100));
            double f1 = Math.log10(errors_points_over_10.get(i) / errors_points_over_10.get(i - 1)) / denom;
            double f2 = Math.log10(errors_points_over_4.get(i) / errors_points_over_4.get(i - 1)) / denom;
            double f3 = Math.log10(errors_points_over_2.get(i) / errors_points_over_2.get(i - 1)) / denom;
            System.out.println(padFloat(f1) + padFloat(f2) + padFloat(f3));
            i++;
        }


    }

    private static double exactSolution(double x) {
        return Math.exp(-(x - 10.) * (x - 10.));
    }

    private static String padFloat(Double number) {
        return padFloat(number, 2, 6);
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
        return padding(number, width - nOfDecimals - 2) + String.format("%." + nOfDecimals + "f", number) + " ";

    }
}
