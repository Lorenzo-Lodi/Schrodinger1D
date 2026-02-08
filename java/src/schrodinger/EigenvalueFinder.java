package schrodinger;

import schrodinger.grid.Grid;
import schrodinger.integrator.Integrator;
import schrodinger.potential.SchrodingerSystem;
import schrodinger.potential.PhysicalPotential;

public class EigenvalueFinder {
    private static final double TARGET_RELATIVE_ERROR = 0.d * Math.ulp(1.d); // change for single-precision float
    private static final int MAXIMUM_NUMBER_OF_BISECTIONS = 60; // reduces error by 2**n
    private final Integrator integrator;
    private final SchrodingerSystem system;

    public EigenvalueFinder(SchrodingerSystem system, Integrator integrator) {
        this.system = system;
        this.integrator = integrator;
    }

    private void propagatePsiLeftToRight(EnergyLevel level) {
        // first (leftmost) point
        level.psi[0] = 0;

        // second point
        level.psi[1] = 0.0000001d; // arbitrary initial value
        system.setEnergy(level.energy);
        for (int n = 1; n < system.getGrid().getNumberOfPoints() - 1; n++) {
            level.psi[n + 1] = integrator.propagateForward(level.psi, n, system);
        }
    }

    /**
     * Locates the energy interval [lowerBound, upperBound] containing the
     * state with 'nOfDesiredNodes' nodes.
     */
    private EnergyLevel findInitialEnergyBracket(int nOfDesiredNodes) {
        Grid grid = system.getGrid();
        EnergyLevel bounds = new EnergyLevel(grid);
        bounds.numberOfNodes = nOfDesiredNodes;

        // 1. Scan the *Effective Potential* U_tilde(y) for the minimum
        int minIndex = 0;
        double uMin = Double.MAX_VALUE;

        // We scan the uniform y-grid
        int nPoints = grid.getNumberOfPoints();
        for (int i = 1; i < nPoints - 1; i++) {
            double y = grid.getYValue(i);
            // Use the effective potential that includes mapping corrections
            double val = system.UTilde(y);
            if (val < uMin) {
                uMin = val;
                minIndex = i;
            }
        }
        bounds.lowerBound = uMin;

        // 2. Estimate Step Size (Energy Scale)
        double step;
        boolean harmonicSuccess = false;

        // Attempt A: Harmonic Curvature of U_tilde on uniform y-grid
        if (minIndex > 0 && minIndex < nPoints - 1) {
            double hy = grid.getStepSizeYCoordinate();
            double yMin = grid.getYValue(minIndex);

            // U_tilde values at minimum and neighbors
            double u0 = system.UTilde(yMin);
            double uL = system.UTilde(yMin - hy);
            double uR = system.UTilde(yMin + hy);

            // Curvature K_y = d^2(U_tilde)/dy^2
            double k_y = (uR - 2.0 * u0 + uL) / (hy * hy);

            // Effective Mass in y-space
            // The kinetic term is - (1 / (2 * m * g^2)) * d^2/dy^2
            // The equation is phi'' = -2m * g^2 * (E - U_tilde) phi
            // So effective mass M_eff = m * g^2(y)
            double g = grid.g(yMin);
            double m_eff = system.getMass() * g * g;

            if (k_y > 1e-15) {
                // Harmonic oscillator: omega = sqrt(K / M)
                step = Math.sqrt(k_y / m_eff);
                harmonicSuccess = true;
            } else {
                step = 0;
            }
        } else {
            step = 0;
        }

        // Attempt B: Particle-in-a-Box (Fallback)
        if (!harmonicSuccess) {
            double L = grid.r(grid.getLastYValue()) - grid.r(grid.getFirstYValue());
            step = (Math.PI * Math.PI) / (2.0 * system.getMass() * L * L);
            if (step < 1e-12) step = 1e-4; // Safety floor
        }

        // 3. Exponential Scan
        double currentEnergy = uMin + step;
        int maxIterations = 100;

        for (int i = 0; i < maxIterations; i++) {
            // Update the energy in the system
            system.setEnergy(currentEnergy);
            bounds.energy = currentEnergy;

            // Propagate using the system (which now has the correct energy)
            propagatePsiLeftToRight(bounds);
            int nodes = bounds.countNumberOfNodes();

            if (nodes > nOfDesiredNodes) {
                bounds.upperBound = currentEnergy;
                return bounds;
            } else {
                bounds.lowerBound = currentEnergy;
                step *= 2.0;
                currentEnergy += step;
            }
        }

        throw new RuntimeException("Failed to bracket energy level.");
    }


    public EnergyLevel findEigenvalueByBisection(int nOfDesiredNodes) {

        EnergyLevel level = this.findInitialEnergyBracket(nOfDesiredNodes);

        // now we can bisect the energy
        for (level.numberOfBisections = 1; level.numberOfBisections <= MAXIMUM_NUMBER_OF_BISECTIONS; level.numberOfBisections++) {
            level.energy = (level.upperBound + level.lowerBound) * 0.5d;
            propagatePsiLeftToRight(level);
            if ((level.upperBound - level.lowerBound) / Math.abs(level.energy) < TARGET_RELATIVE_ERROR) {
                break;
            }
            if (level.countNumberOfNodes() > nOfDesiredNodes) {
                level.upperBound = level.energy;
            } else {
                level.lowerBound = level.energy;
            }
        }
        system.setEnergy(level.energy);
        level.normalizePsi();
        level.perturbativeCorrectionToEnergy = integrator.computePerturbativeCorrection(level, system);

        return level;
    }

    // TODO implement findEigenvalueBySecant


}
