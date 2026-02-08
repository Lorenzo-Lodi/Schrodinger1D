package schrodinger.grid;

public class MappingSurkus implements MappingStrategy {
    private final double rRef;
    private final double alpha;

    public MappingSurkus(double rRef, double alpha) {
        this.rRef = rRef;
        this.alpha = alpha;
    }

    @Override
    public double y(double r) {
        return (Math.pow(r / rRef, alpha) - 1.0d) / (Math.pow(r / rRef, alpha) + 1.0d);
    }

    @Override
    public double r(double y) {
        return rRef * Math.pow((1 + y) / (1 - y), 1.0d / alpha);
    }

    @Override
    public double g(double y) {
        return (2.0d * rRef * Math.pow(1 + y, -1.0d + 1.0d / alpha))
                / (alpha * Math.pow(1.0d - y, 1.0d + 1.0d / alpha));
    }

    @Override
    public double F(double y) {
        return (1.0d - 1.0d / Math.pow(alpha, 2)) / Math.pow(1 - y * y, 2);
    }

}
