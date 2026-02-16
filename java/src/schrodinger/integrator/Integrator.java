package schrodinger.integrator;

import schrodinger.QuantumState;
import schrodinger.potential.SchrodingerSystem;

public interface Integrator {

    double propagate(double[] psi, int n, QuantumState state, Direction direction);

    // TODO should probably go into a separate interface.
    default void computePerturbativeCorrection(QuantumState state) {
    }

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
