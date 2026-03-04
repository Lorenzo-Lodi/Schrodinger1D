package schrodinger.potential;

public class PhysicalPotentialPolynomial implements PhysicalPotential {
    private final double r0;
    private final double[] coefficients;

    public PhysicalPotentialPolynomial(double r0, double[] coefficients) {
        this.r0 = r0;
        this.coefficients = coefficients;
    }

    public double value(double r) {
        double ret = 0;
        for (int i = 0; i < coefficients.length; i++) {
            ret = ret + coefficients[i] * Math.pow(r - r0, i);
        }
        return ret;
    }
}
