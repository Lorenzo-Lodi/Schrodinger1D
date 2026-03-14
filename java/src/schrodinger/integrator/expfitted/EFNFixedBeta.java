package schrodinger.integrator.expfitted;

// Exponentially-Fitted Numerov Method with Fixed β
public class EFNFixedBeta extends ExponentiallyFittedAbstract {

    public double getBeta(double Z) {
        return 1.0 / 12.0;
    }

    public double getGamma(double Z, double beta) {

        if (Math.abs(Z) < 0.05) {
            // Extended Taylor series for gamma(Z) around Z = 0
            // Valid for both Z > 0 (oscillatory) and Z < 0 (exponential)
            double Z2 = Z * Z;
            double Z3 = Z2 * Z;
            double Z4 = Z2 * Z2;
            return (5.0 / 6.0) - (Z2 / 240.0) + (11.0 * Z3 / 60480.0) - (13.0 * Z4 / 3628800.0);
        }

        double sqrtModZ = Math.sqrt(Math.abs(Z));
        double c;

        if (Z > 0) { // Oscillatory Region (Q > 0)
            c = Math.cos(sqrtModZ);
        } else { // Exponential Region (Q < 0)
            c = Math.cosh(sqrtModZ);
        }
        return 2.0 * (1.0 - c * (1.0 + beta * Z)) / Z;
    }
}


