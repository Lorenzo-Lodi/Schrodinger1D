package schrodinger;

import schrodinger.grid.Grid;
import schrodinger.integrator.Integrator;
import schrodinger.potential.TransformedQFunction;
import schrodinger.potential.PhysicalPotential;

public class EigenvalueFinder {
    private static final double TARGET_RELATIVE_ERROR = 0.d * Math.ulp(1.d); // change for single-precision float
    private static final int MAXIMUM_NUMBER_OF_BISECTIONS = 60; // reduces error by 2**n
    private final Grid grid;
    private final PhysicalPotential physicalPotential;
    private final double mass;
    private final Integrator integrator;

    public EigenvalueFinder(Grid grid, PhysicalPotential physicalPotential, double mass, Integrator integrator) {
        this.grid = grid;
        this.physicalPotential = physicalPotential;
        this.mass = mass;
        this.integrator = integrator;
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

        TransformedQFunction transformedQFunction = new TransformedQFunction(physicalPotential, mass, level.energy, grid);
        level.normalizePsi();
        level.perturbativeCorrectionToEnergy = integrator.computePerturbativeCorrection(level, transformedQFunction);

        return level;
    }

    // TODO implement findEigenvalueBySecant

    private void propagatePsiLeftToRight(EnergyLevel level) {
        // first (leftmost) point
        level.psi[0] = 0;

        // second point
        level.psi[1] = 0.0000001d; // arbitrary initial value

        TransformedQFunction transformedQFunction = new TransformedQFunction(physicalPotential, mass, level.energy, grid);
        for (int n = 1; n < grid.getNumberOfPoints() - 1; n++) {
            level.psi[n + 1] = integrator.propagateForward(level.psi, n, transformedQFunction);
        }
    }

    /**
     * Locates the energy interval [lowerBound, upperBound] that contains the
     * eigenvalue with the specified number of nodes.
     * <p>
     * This method replaces hard-coded "magic numbers" with a physics-based
     * heuristic. It estimates the energy scale using the harmonic curvature
     * of the potential at its minimum or the confinement scale of the grid.
     *
     * @param nOfDesiredNodes The target number of nodes for the eigenstate.
     * @return An EnergyLevel object with lowerBound and upperBound populated.
     */
    private EnergyLevel findInitialEnergyBracket(int nOfDesiredNodes) {
        EnergyLevel bounds = new EnergyLevel(grid);
        bounds.numberOfNodes = nOfDesiredNodes;

        // 1. Scan the grid to find the global minimum of the potential
        int minIndex = 0;
        double vMin = Double.MAX_VALUE;

        // We scan the physical potential V(r) to find the absolute energy floor
        for (int i = 0; i < grid.getNumberOfPoints(); i++) {
            double r = grid.r(grid.getYValue(i));
            double v = physicalPotential.value(r);
            if (v < vMin) {
                vMin = v;
                minIndex = i;
            }
        }
        bounds.lowerBound = vMin;

        // 2. Estimate the initial search step (Energy Scale)
        double step;
        boolean harmonicSuccess = false;

        // Attempt A: Harmonic Approximation (Local)
        // We estimate curvature (k) using the uniform computational y-grid
        if (minIndex > 0 && minIndex < grid.getNumberOfPoints() - 1) {
            double hy = grid.getStepSizeYCoordinate();

            // Potential values at the minimum and neighbors
            double v0 = physicalPotential.value(grid.r(grid.getYValue(minIndex)));
            double vL = physicalPotential.value(grid.r(grid.getYValue(minIndex - 1)));
            double vR = physicalPotential.value(grid.r(grid.getYValue(minIndex + 1)));

            // Second derivative d2V/dy2 on uniform grid
            double vyy = (vR - 2.0 * v0 + vL) / (hy * hy);

            // Map the curvature from y-space back to r-space: k_r = Vyy / (dr/dy)^2
            // dr/dy (central difference)
            double rL = grid.r(grid.getYValue(minIndex - 1));
            double rR = grid.r(grid.getYValue(minIndex + 1));
            double drdy = (rR - rL) / (2.0 * hy);

            double k = vyy / (drdy * drdy);

            // If k > 0, we have a stable well; estimate step as h_bar * omega
            if (k > 1e-12) {
                step = Math.sqrt(k / mass); // Assuming units where h_bar = 1
                harmonicSuccess = true;
            } else {
                step = 0; // Fallback
            }
        } else {
            step = 0; // Fallback
        }

        // Attempt B: Particle-in-a-Box (Global Fallback)
        // Used if the potential is flat (k=0) or inverted (k<0) at the minimum
        if (!harmonicSuccess) {
            double L = grid.r(grid.getLastYValue()) - grid.r(grid.getFirstYValue());
            // E_ground_state = pi^2 / (2 * mass * L^2)
            step = (Math.PI * Math.PI) / (2.0 * mass * L * L);

            // Absolute safety floor to avoid zero step
            if (step < 1e-12) step = 1e-4;
        }

        // 3. Exponential Scan: Double the step until we bracket the state
        // currentEnergy starts one step above the minimum
        double currentEnergy = vMin + step;
        int maxIterations = 100; // Sufficient to cover ~30 orders of magnitude

        for (int i = 0; i < maxIterations; i++) {
            bounds.energy = currentEnergy;

            // propagatePsiLeftToRight updates bounds.psi and node count
            propagatePsiLeftToRight(bounds);
            int nodes = bounds.countNumberOfNodes();

            if (nodes > nOfDesiredNodes) {
                // Success: This energy has too many nodes, so the bracket is closed
                bounds.upperBound = currentEnergy;
                return bounds;
            } else {
                // Still too low: update lowerBound and accelerate upwards
                bounds.lowerBound = currentEnergy;
                step *= 2.0;
                currentEnergy += step;
            }
        }

        throw new RuntimeException("Failed to bracket energy level. Potential may be unbound or mass is too small.");
    }


}
