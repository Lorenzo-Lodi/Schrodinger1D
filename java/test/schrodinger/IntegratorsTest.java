package schrodinger;

import org.junit.jupiter.api.Test;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.*;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialHarmonic;
import schrodinger.potential.SchrodingerSystem;


public class IntegratorsTest {

    @Test
    public void propagate_test() {
        PhysicalPotential potential = new PhysicalPotentialHarmonic(1.2, 1.1);
        int nOfPoints = 10;

        Grid grid = GridFactory.generateUniformGrid(1.2, 1.6, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, 1.5, grid);
        QuantumState state = new QuantumState(system);
        state.energy = -1000.;
        state.psi[0] = 1.0;
        state.psi[1] = 1.5;

        System.out.println(new TaylorThreePoints().propagate(state.psi, 1, state, Integrator.Direction.FORWARD));
        System.out.println(new Numerov().propagate(state.psi, 1, state, Integrator.Direction.FORWARD));
        System.out.println(new Vignoli().propagate(state.psi, 1, state, Integrator.Direction.FORWARD));
        System.out.println(new ExponentiallyFitted().propagate(state.psi, 1, state, Integrator.Direction.FORWARD));
    }

    @Test
    public void efn_test() {

        ExponentiallyFittedAbstract int1 = new Vignoli();
        ExponentiallyFittedAbstract int2 = new ExponentiallyFitted();

        for (int k = 0; k <= 1; k++) {
            for (int i = -2; i <= 12; i++) {
                int exp;
                if (k == 0) {
                    exp = -i;
                } else {
                    exp = -10 + i;
                }
                double z = Math.pow(10., exp);
                if (k == 0) z = -z;
                System.out.println(padFloat(z, 15, 20) + " " + padFloat(int1.getBeta(z)) + " " + padFloat(int2.getBeta(z))
                        + "          "
                        + padFloat(int1.getGamma(z, int1.getBeta(z))) + " " + padFloat(int2.getGamma(z, int2.getBeta(z))));
            }
        }

    }


    @Test
    public void efn_test2() {

        ExponentiallyFittedAbstract int1 = new Vignoli();
        ExponentiallyFittedAbstract int2 = new ExponentiallyFitted();

        for (int k = -1000; k <= 1000; k++) {
            double z = ((double) k) / 10.;
            System.out.println(padFloat(z, 15, 20) + " " + padFloat(int1.getBeta(z)) + " " + padFloat(int2.getBeta(z))
                    + "          "
                    + padFloat(int1.getGamma(z, int1.getBeta(z))) + " " + padFloat(int2.getGamma(z, int2.getBeta(z))));
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
