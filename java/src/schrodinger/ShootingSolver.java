package schrodinger;

import schrodinger.grid.Grid;
import schrodinger.integrator.Integrator;
import schrodinger.potential.SchrodingerSystem;

public class ShootingSolver {
    private static final double TARGET_ABSOLUTE_ERROR = 1e-8;
    private static final int MAXIMUM_NUMBER_OF_BISECTIONS = 50; // reduces error by 2**n
    private final Integrator integrator;
    private final SchrodingerSystem system;

    public ShootingSolver(SchrodingerSystem system, Integrator integrator) {
        this.system = system;
        this.integrator = integrator;
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
        double energyScale;
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

            if (k_y > 1e-15) {
                // Harmonic oscillator: omega = sqrt(K / M)
                energyScale = Math.sqrt(k_y / system.getMass());
                harmonicSuccess = true;
            } else {
                energyScale = 0;
            }
        } else {
            energyScale = 0;
        }

        // Attempt B: Particle-in-a-Box (Fallback)
        if (!harmonicSuccess) {
            double L = grid.getLastYValue() - grid.getFirstYValue();
            energyScale = (Math.PI * Math.PI) / (2.0 * system.getMass() * L * L);
            if (energyScale < 1e-12) energyScale = 1e-4; // Safety floor
        }

        // 3. Exponential Scan to find upper bound to the energy
        double currentEnergy = uMin + energyScale;
        int maxIterations = 100;

        for (int i = 0; i < maxIterations; i++) {
            // Update the energy in the system
            system.setEnergy(currentEnergy);
            bounds.energy = currentEnergy;
            int nodes = countNodesRobust(bounds);

            if (nodes > nOfDesiredNodes) {
                bounds.upperBound = currentEnergy;
                return bounds;
            } else {
                bounds.lowerBound = currentEnergy;
                energyScale *= 2.0;
                currentEnergy += energyScale;
            }
        }

        throw new RuntimeException("Failed to bracket energy level.");
    }


    /**
     * Counts the number of nodes for a given energy using a stable
     * bidirectional shooting method.
     *
     * Strategy:
     * 1. Find the classical turning point (matching point).
     * 2. Shoot Forward from Left -> Match.
     * 3. Shoot Backward from Right -> Match.
     * 4. Match signs at the junction.
     * 5. Sum nodes.
     */
    /**
     * Counts nodes using bidirectional shooting (stable for large grids).
     */
    private int countNodesRobust(EnergyLevel level) {
        int nPoints = system.getGrid().getNumberOfPoints();
        int matchIndex = findMatchingIndex(level.energy);

        // --- Shoot Forward: Left boundary -> Match ---
        level.psi[0] = 0.0;
        level.psi[1] = 1e-16;

        int nodesLeft = 0;
        for (int n = 1; n < matchIndex; n++) {
            level.psi[n + 1] = integrator.propagate(level.psi, n, system, Integrator.Direction.FORWARD);
            if (level.psi[n] * level.psi[n + 1] < 0.0) {
                nodesLeft++;
            }
        }

        // --- Shoot Backward: Right boundary -> Match ---
        level.psi[nPoints - 1] = 0.0;
        level.psi[nPoints - 2] = 1e-16;

        int nodesRight = 0;
        for (int n = nPoints - 2; n > matchIndex; n--) {
            level.psi[n - 1] = integrator.propagate(level.psi, n, system, Integrator.Direction.BACKWARD);
            if (level.psi[n] * level.psi[n - 1] < 0.0) {
                nodesRight++;
            }
        }

        // Total nodes: simply sum from both segments
        return nodesLeft + nodesRight;
    }

    private int findMatchingIndex(double energy) {
        Grid grid = system.getGrid();
        // Scan from right to left, within safe range
        for (int i = grid.getNumberOfPoints() - 3; i >= 2; i--) {
            if (system.UTilde(grid.getYValue(i)) <= energy) {
                return i;
            }
        }

        // Fallback to midpoint if always classically forbidden
        return grid.getNumberOfPoints() / 2;
    }


    public EnergyLevel findEigenvalueByBisection(int nOfDesiredNodes) {

        EnergyLevel level = this.findInitialEnergyBracket(nOfDesiredNodes);

        // now we can bisect the energy
        for (level.numberOfBisections = 1; level.numberOfBisections <= MAXIMUM_NUMBER_OF_BISECTIONS; level.numberOfBisections++) {
            double mid = (level.lowerBound + level.upperBound) * 0.5;
            system.setEnergy(mid);
            level.energy = mid;
            int nodes = countNodesRobust(level);

            if (Math.abs((level.upperBound - level.lowerBound)) < TARGET_ABSOLUTE_ERROR) {
                break;
            }
            if (nodes > nOfDesiredNodes) {
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
