package schrodinger.potential;

import schrodinger.grid.Grid;

public class SchrodingerSystem {
    private final PhysicalPotential physicalPotential;
    private final double mass;
    private final Grid grid;
    private final double qMin;
    private final static double hCritical = 2.5;
    private boolean isCapPotential;

    public SchrodingerSystem(PhysicalPotential physicalPotential, double mass, Grid grid) {
        this.physicalPotential = physicalPotential;
        this.mass = mass;
        this.grid = grid;
        qMin = -Math.pow(hCritical / grid.getStepSizeYCoordinate(), 2);
    }

    /**
     * Potential U(r) for the equation: -(hbar^2 / 2m) ψ''(r) + U(r) ψ(r) = E ψ(r)
     */
    public double U(double r) {
        return physicalPotential.value(r);
    }

    private double UTilde(double y) {
        double u = physicalPotential.value(grid.r(y));
        double gy = grid.g(y);
        double factor = grid.F(y) / (gy * gy);
        return u - factor / (2.d * mass);
    }

    public double UTildeAtGridPoint(double i) {
        double y = grid.yAtGridPoint(i);
        return UTilde(y);
    }

    public double getMass() {
        return mass;
    }

    public Grid getGrid() {
        return grid;
    }

    public double getQMin() {
        return qMin;
    }

    public boolean isCapPotential() {
        return isCapPotential;
    }

    public void setCapPotential(boolean capPotential) {
        isCapPotential = capPotential;
    }

}
