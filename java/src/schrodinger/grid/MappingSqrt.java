package schrodinger.grid;

public class MappingSqrt implements MappingStrategy {
    private final double rRef;

    public MappingSqrt(double rRef) {
        this.rRef = rRef;
    }

    @Override
    public double y(double r) {
        return 2.d * rRef * (Math.sqrt(1 + r / rRef) - 1.d);
    }

    @Override
    public double r(double y) {
        return y * (1. + y / (4.0d * rRef));
    }

    @Override
    public double g(double y) {
        return 1.0d + y / (2.0d * rRef);
    }

    @Override
    public double F(double y) {
        return -3.0 / Math.pow(8.0d * rRef + 4.0d * y, 2);
    }

}
