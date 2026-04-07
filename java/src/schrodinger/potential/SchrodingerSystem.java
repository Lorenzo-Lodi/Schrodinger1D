package schrodinger.potential;

import schrodinger.grid.Grid;

public class SchrodingerSystem {
    private final PhysicalPotential physicalPotential;
    private final double mass;
    private final Grid grid;
    private final double qMin;

    // For now we use a unique hCritical for all integrators. In reality some integrators are very sensitive
    // (PC6 is the most sensitive and requires hCritical=2.5) and some much less (Obrechkoff6 requires hCritical=5.5)
    // That said, this corresponds to a factor (5.5/2.5)^2 = 4.84 in the capping value, which is not huge.
    // These are the critical values (tested with Lennard-Jones potential) with 0.5 accuracy (only for integrators
    // which need capping):
//    Obrechkoff6            5.5
//    Stormer5               3.5
//    Stormer6               3.5
//    Stormer8i              3.5
//    PredictorCorrector8i2  3.5
//    Numerov                3.0
//    EFNFixedBeta           3.0
//    PredictorCorrector8i1  3.0
//    PredictorCorrector6    2.5
    private final static double hCritical = 2.5;

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

}
