package schrodinger.integrator;

public class PadeNumerov4 extends ExponentiallyFittedAbstract {

    @Override
    public double getBeta(double Z) {
        return 1.0 / 12.0 + Z / 240.0 + Z * Z / 6048.0;
    }

    @Override
    public double getGamma(double Z, double beta) {
        return 1.0 - 2.0 * beta;
    }
}
