package schrodinger.solver;

import schrodinger.QuantumState;
import schrodinger.grid.Grid;
import schrodinger.integrator.Integrator;
import schrodinger.potential.SchrodingerSystem;

public class ShootingSolver {
    private static final double TARGET_ABSOLUTE_ERROR = 1e-13;
    private static final int MAXIMUM_NUMBER_OF_BISECTIONS = 60; // reduces error by 2**n
    private static final double PSI_MAX = 1e30; // stop integrating forward if wave function exceeds this value
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
    private QuantumState findInitialEnergyBracket(int nOfDesiredNodes) {
        QuantumState level = new QuantumState(system);
        level.numberOfNodes = nOfDesiredNodes;

        double energyScale = findSystemEnergyScale(level);

        // 3. Exponential Scan to find upper bound to the energy
        double currentEnergy = level.energy;

        for (int i = 0; i < MAXIMUM_NUMBER_OF_BISECTIONS; i++) {
            // Update the energy in the system
            level.energy = currentEnergy;
            int nodes = countNodes(level);

            if (nodes > nOfDesiredNodes) {
                level.upperBound = currentEnergy;
                level.numberOfNodesUpperBound = nodes;
                return level;
            } else {
                level.lowerBound = currentEnergy;
                level.numberOfNodesLowerBound = nodes;
                energyScale *= 2.0;
                currentEnergy += energyScale;
            }
        }

        throw new RuntimeException("Failed to bracket energy level.");
    }

    private Double findSystemEnergyScale(QuantumState level) {
        Grid grid = system.getGrid();

        // 1. Scan the *Effective Potential* U_tilde(y) for the minimum
        int minIndex = 0;
        double uMin = Double.MAX_VALUE;

        // We scan the uniform y-grid
        int nPoints = grid.getNumberOfPoints();
        for (int i = 1; i < nPoints - 1; i++) {
            // Use the effective potential that includes mapping corrections
            double val = system.UTildeValueAt(i);
            if (val < uMin) {
                uMin = val;
                minIndex = i;
            }
        }
        level.lowerBound = uMin;
        level.numberOfNodesLowerBound = 0; // Should be always correct

        // 2. Estimate Step Size (Energy Scale)
        Double energyScale = null;

        // Attempt A: Harmonic Curvature of U_tilde on uniform y-grid
        if (minIndex > 0 && minIndex < nPoints - 1) {

            // U_tilde values at minimum and neighbors
            double u0 = system.UTildeValueAt(minIndex);
            double uL = system.UTildeValueAt(minIndex - 1);
            double uR = system.UTildeValueAt(minIndex + 1);

            // second derivative  d^2(U_tilde)/dy^2
            double hy = grid.getStepSizeYCoordinate();
            double der2 = (uR - 2.0 * u0 + uL) / (hy * hy);

            if (der2 > 1e-15) {
                // Harmonic oscillator: omega = sqrt(K / M)
                energyScale = Math.sqrt(der2 / system.getMass());
            }
        }

        // Attempt B: Particle-in-a-Box (Fallback)
        if (energyScale == null) {
            double L = grid.getLastYValue() - grid.getFirstYValue();
            energyScale = (Math.PI * Math.PI) / (2.0 * system.getMass() * L * L);
        }

        level.energy = uMin + energyScale;

        return energyScale;

    }

    public QuantumState findEigenvalueByBisection(int nOfDesiredNodes) {

        QuantumState level = this.findInitialEnergyBracket(nOfDesiredNodes);
        refineByBisection(level, nOfDesiredNodes, TARGET_ABSOLUTE_ERROR);
        integrator.computePerturbativeCorrection(level);

        return level;
    }

    private void refineByBisection(QuantumState level, int nOfDesiredNodes, double targetAbsoluteError) {

        for (level.numberOfBisections = 1; level.numberOfBisections <= MAXIMUM_NUMBER_OF_BISECTIONS; level.numberOfBisections++) {
            double mid = (level.lowerBound + level.upperBound) * 0.5;
            level.energy = mid;
            int nodes = countNodes(level);

            if (nodes > nOfDesiredNodes) {
                level.upperBound = level.energy;
                level.numberOfNodesUpperBound = nodes;
            } else {
                level.lowerBound = level.energy;
                level.numberOfNodesLowerBound = nodes;
            }

            if (Math.abs((level.upperBound - level.lowerBound)) < targetAbsoluteError &&
                    level.numberOfNodesLowerBound == nOfDesiredNodes
                    && level.numberOfNodesUpperBound == nOfDesiredNodes + 1) {
                break;
            }

        }

    }


    // TODO WIP on bi-directional matching method
    public QuantumState findEigenvalueHybridMethod(int nOfDesiredNodes) {
        QuantumState level = this.findInitialEnergyBracket(nOfDesiredNodes);
        refineByBisection(level, nOfDesiredNodes, 1e-3);



        integrator.computePerturbativeCorrection(level);

        return level;
    }


    // TODO implement findEigenvalueBySecant

    private int countNodes(QuantumState level) {
        int nPoints = system.getGrid().getNumberOfPoints();

        // --- Shoot Forward
        level.psi[0] = 0.0;
        level.psi[1] = 1e-16;

        int nodes = 0;
        for (int n = 1; n < nPoints - 1; n++) {
            level.psi[n + 1] = integrator.propagate(level.psi, n, level, Integrator.Direction.FORWARD);
            if (level.psi[n] * level.psi[n + 1] < 0.0) {
                nodes++;
            }

            // If we are deep in the classically-forbidded region and the wavefunction is blowing up, we can stop early
            if (level.QTildeValueAt(n) > 2.0 * level.energy && Math.abs(level.psi[n + 1]) > PSI_MAX) {
                break;
            }
        }
        return nodes;
    }

    private int findMatchingIndex(double energy) {
        Grid grid = system.getGrid();
        // Scan from right to left until we reach the classically-allowed region
        for (int i = grid.getNumberOfPoints() - 3; i >= 2; i--) {
            if (system.UTildeValueAt(i) <= energy) {
                return i;
            }
        }

        // Fallback to midpoint if always classically forbidden
        return grid.getNumberOfPoints() / 2;
    }

}
