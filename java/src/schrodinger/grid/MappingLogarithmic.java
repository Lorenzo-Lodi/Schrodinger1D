package schrodinger.grid;

public class MappingLogarithmic implements MappingStrategy {
    double rRef;

    public MappingLogarithmic(double rRef) {
        this.rRef = rRef;
    }

    @Override
    public double y(double r) {
        return rRef * Math.log(1 + r / rRef);
    }

    @Override
    public double r(double y) {
        return rRef * (-1. + Math.exp(y / rRef));
    }

    @Override
    public double g(double y) {
        return Math.exp(y / rRef);
    }

    @Override
    public double F(double y) {
        return -1. / (4. * rRef * rRef);
    }

}
