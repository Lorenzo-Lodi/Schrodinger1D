package schrodinger.potential;

import schrodinger.grid.Grid;

public class SchrodingerSystem {
    private final PhysicalPotential physicalPotential;
    private final double mass;
    private double energy;
    private final Grid grid;

    public SchrodingerSystem(PhysicalPotential physicalPotential, double mass, double energy, Grid grid) {
        this.physicalPotential = physicalPotential;
        this.mass = mass;
        this.energy = energy;
        this.grid = grid;
    }

    public double U(double r) {
        return physicalPotential.value(r);
    }

    public double UTilde(double y) {
        double u = physicalPotential.value(grid.r(y));
        double factor = grid.F(y) / Math.pow(grid.g(y), 2);
        return u - factor / (2.d * mass);
    }

    public double Q(double r) {
        return 2.d * mass * (energy - physicalPotential.value(r));
    }

    /**
     * Transformed Q-function Q̃(y) for the mapped equation: ϕ''(y) = -Q̃(y)ϕ(y)
     * MSL Eq. (8): Q̃(y) = g²(y)·Q(r(y)) + F(y)
     */
    public double QTilde(double y) {
        return Math.pow(grid.g(y), 2) * Q(grid.r(y)) + grid.F(y);
    }

    public double getMass() {
        return mass;
    }

    public Grid getGrid() {
        return grid;
    }

    public double getEnergy() {
        return energy;
    }

    public void setEnergy(double energy) {
        this.energy = energy;
    }

}
