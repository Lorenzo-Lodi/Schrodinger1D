package schrodinger.integrator;

public class Vignoli extends ExponentiallyFittedAbstract {

    double getBeta(double Z) {
        return 1.0 / 12.0;
    }

    double getGamma(double Z, double beta) {

        if (Math.abs(Z) < 1e-4) {
            return 5.0 / 6.0; // Standard Numerov value
        }

        double sqrtModZ = Math.sqrt(Math.abs(Z));
        double c;

        if (Z > 0) { // Oscillatory Region (Q > 0)
            c = 1. - Math.cos(sqrtModZ);
        } else { // Exponential Region (Q < 0)
            c = 1. - Math.cosh(sqrtModZ);
        }
        return (2.0 - 2.0 * (1. - c) * (1.0 + beta * Z)) / Z;
    }
}


