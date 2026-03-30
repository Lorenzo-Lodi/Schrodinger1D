package schrodinger.potential;

public class PhysicalPotentialLennardJones implements PhysicalPotential {
    private final double rmin;
    private final double wellDepth;
    private final int exponent;
    private final double sigma;

    public PhysicalPotentialLennardJones(double rmin, double wellDepth, int exponent) {
        this.rmin = rmin;
        this.wellDepth = wellDepth;
        this.exponent = exponent;
        this.sigma = rmin * Math.pow(2., -1. / exponent);
    }

    public double value(double r) {
        double x = Math.pow(sigma / r, exponent);
        return 4. * wellDepth * x * (x - 1.) + wellDepth;
    }

}
