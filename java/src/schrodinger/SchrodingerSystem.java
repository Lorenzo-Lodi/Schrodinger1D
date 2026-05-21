package schrodinger;

import schrodinger.grid.Grid;
import schrodinger.potential.PhysicalPotential;

import static schrodinger.PhysicalConstants.toInverseCm;

public class SchrodingerSystem {
    private final PhysicalPotential physicalPotential;
    private final double mass;
    private final Grid grid;
    private final double qMin;
    public double UTildeMinimumGridIndex;
    public double UTildeMinimumGridValue;
    private Double energyScale = null;
    private FractionalGridCache cacheUTilde;

    // For now we use a unique hCritical for all integrators. In reality some integrators are very sensitive
    // (PC6 is the most sensitive and requires hCritical=2.5) and some much less (Obrechkoff6 requires hCritical=5.5)
    // That said, this corresponds to a factor (5.5/2.5)^2 = 4.84 in the capping value, which is not huge.
    // The vast majority require 3.5, so the factor (3.5/2.5)^2 = 1.96, ie setting to 2.5 we cap to HALF the value we
    // could be capping.
    // These are the critical values (tested with Lennard-Jones potential) with 0.5 accuracy (only for integrators
    // which need capping):
    //    Obrechkoff6            5.5
    //    Cowell5                3.5
    //    Cowell6                3.5
    //    Cowell8                3.5
    //    PC8i2                  3.5
    //    Numerov                3.0
    //    EFNFixedBeta           3.0
    //    PC8i1                  3.0
    //    PC6                    2.5
    private final static double hCriticalForbidden = 2.5;
    public double hCriticalAllowed = 0;

    public SchrodingerSystem(PhysicalPotential physicalPotential, double mass, Grid grid) {
        this.physicalPotential = physicalPotential;
        this.mass = mass;
        this.grid = grid;
        qMin = -Math.pow(hCriticalForbidden / grid.getStepSizeYCoordinate(), 2);
        this.estimateEnergyScale(); // Pre-compute the energy scale
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

    public double Ucapped() {
        return -qMin / (2. * mass);
    }


    public double estimateEnergyScale() {
        if (energyScale != null) {
            OutputManager.write(String.format("The energy scale was set to %23.14f (%25.6f cm-1) ",
                    energyScale, toInverseCm(energyScale)));
            return energyScale;
        }

        // 1. Scan the *Effective Potential* U_tilde(y) for the minimum
        int minIndex = 0;
        double uMin = Double.MAX_VALUE;

        // We scan the uniform y-grid
        int nPoints = grid.getNumberOfPoints();
        for (int i = 1; i < nPoints - 1; i++) {
            // Use the effective potential that includes mapping corrections
            double val = UTildeAtGridPoint(i);
            if (val < uMin) {
                uMin = val;
                minIndex = i;
            }
        }
        UTildeMinimumGridIndex = minIndex;
        UTildeMinimumGridValue = uMin;

        OutputManager.write(String.format("I scanned the potential and found a minimum value %23.14f (%25.6f cm-1) for i = %d",
                uMin, toInverseCm(uMin), minIndex));

        double uMaxRight = -(0.5 * Double.MAX_VALUE);
        for (int i = minIndex + 1; i < nPoints - 1; i++) {
            double val = UTildeAtGridPoint(i);
            if (val > uMaxRight) {
                uMaxRight = val;
            }
        }
        OutputManager.write(String.format("I scanned the potential and found a maximum value %23.14f (%25.6f cm-1) right of the minimum",
                uMaxRight, toInverseCm(uMaxRight)));
        double kMax = Math.sqrt(Math.max(0, 2. * mass * (uMaxRight - uMin)));
        // Minimum step size for all states up to uMaxRight (for Morse-like potential, the dissociation energy of the potential)
        // For potentials which are not-Morse like, this is wrong.
        hCriticalAllowed = Math.PI / (4. * kMax); // The number at the denominator is somewhat arbitrary, values for 2 to 4 are reasonable.

        // 2. Estimate Step Size (Energy Scale)

        // Attempt A: Harmonic Curvature of U_tilde on uniform y-grid
        String energyScaleMethod = "";
        if (minIndex > 0 && minIndex < nPoints - 1) {

            // U_tilde values at minimum and neighbors
            double u0 = UTildeAtGridPoint(minIndex);
            double uL = UTildeAtGridPoint(minIndex - 1);
            double uR = UTildeAtGridPoint(minIndex + 1);

            // second derivative  d^2(U_tilde)/dy^2
            double hy = grid.getStepSizeYCoordinate();
            double der2 = (uR - 2.0 * u0 + uL) / (hy * hy);

            if (der2 > 1e-15) {
                // Harmonic oscillator: omega = sqrt(K / M)
                energyScale = Math.sqrt(der2 / getMass());
                energyScaleMethod = "harmonic constant at equilibrium";
            }
        }

        // Attempt B: Particle-in-a-Box (Fallback)
        if (energyScale == null) {
            double L = grid.getLastYValue() - grid.getFirstYValue();
            energyScale = (Math.PI * Math.PI) / (2.0 * getMass() * L * L);
            energyScaleMethod = "Particle-in-a-Box";
        }

        OutputManager.write(String.format("The energy scale was set to %23.14f (%25.6f cm-1) using as method: %s",
                energyScale, toInverseCm(energyScale), energyScaleMethod));

        return energyScale;
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

    public double findLeftmostInversionGridpointY(double energy) {
        for (int i = 0; i < grid.getNumberOfPoints(); i++) {
            if (UTildeAtGridPoint(i) <= energy) {
                return grid.yAtGridPoint(i);
            }
        }
        return grid.getFirstYValue();
    }

    public double findRightmostInversionGridpointY(double energy) {
        for (int i = grid.getNumberOfPoints() - 1; i > 0; i--) {
            if (UTildeAtGridPoint(i) <= energy) {
                return grid.yAtGridPoint(i);
            }
        }
        return grid.getLastYValue();
    }

    public void initializeCache(double[] offsets) {
        this.cacheUTilde = new FractionalGridCache(this.grid.getNumberOfPoints(), offsets, this::UTildeAtGridPoint);
    }

}
