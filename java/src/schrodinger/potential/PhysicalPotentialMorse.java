package schrodinger.potential;

public class PhysicalPotentialMorse implements PhysicalPotential {
    private final double rmin;
    private final double a;
    private final double De;

    public PhysicalPotentialMorse(double rmin, double a, double De) {
        this.rmin = rmin;
        this.De = De;
        this.a = a;
    }

    public double value(double r) {
        double x = 1. - Math.exp(-a * (r - rmin));
        return De * x * x;
    }

}
