package schrodinger.integrator;

import java.util.function.DoubleUnaryOperator;

/**
 * Interface for numerical integration methods used to propagate the wavefunction.
 */
public interface Integrator {

    /**
     * Propagates the wavefunction from index n to index n+1 or n-1 depending on direction.
     *
     * @param psi       The wavefunction array containing partially filled values
     * @param n         The current index (point n is assumed to be computed)
     * @param step      The step size in the y-coordinate
     * @param qTilde    Function that returns the Q-tilde value at a given index
     * @param direction The integration direction (FORWARD or BACKWARD)
     * @return The value of the wavefunction at the next point (n+1 or n-1)
     */

    double propagate(double[] psi, double[] currentPsiPrime, int n, double step,
                     DoubleUnaryOperator qTilde,
                     DoubleUnaryOperator qTildePrime,
                     DoubleUnaryOperator qTildeDoublePrime, Direction direction);

    /**
     * Returns the minimum number of previously-computed psi values this integrator
     * needs to propagate one step. A value of 2 means psi[n] and psi[n-1] must both
     * be known; 3 means psi[n], psi[n-1], psi[n-2] must be known; etc.
     *
     * <p>The actual psi array may be longer — this is only the minimum requirement.
     */
    int minHistoryLength();

    /**
     * Returns the global convergence order of this integration method.
     * For a method of order p, the global error scales as O(h^p) as the step size h → 0.
     *
     * <p>This is informative metadata about the algorithm's mathematical properties.
     * It is also used by the test framework to calibrate convergence tolerance checks.
     */
    int globalConvergenceOrder();

    enum Direction {
        FORWARD(1), BACKWARD(-1);

        private final int value;

        Direction(int value) {
            this.value = value;
        }

        public int getValue() {
            return value;
        }
    }

}
