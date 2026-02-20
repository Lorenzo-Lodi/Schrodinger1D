package schrodinger.integrator;

import java.util.function.IntToDoubleFunction;

/**
 * Interface for numerical integration methods used to propagate the wavefunction.
 */
public interface Integrator {

    /**
     * Propagates the wavefunction from index n to index n+1 or n-1 depending on direction.
     * 
     * @param psi The wavefunction array containing partially filled values
     * @param n The current index (point n is assumed to be computed)
     * @param step The step size in the y-coordinate
     * @param qTildeFunction Function that returns the Q-tilde value at a given index
     * @param direction The integration direction (FORWARD or BACKWARD)
     * @return The value of the wavefunction at the next point (n+1 or n-1)
     */
    double propagate(double[] psi, int n, double step, IntToDoubleFunction qTildeFunction, Direction direction);

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
