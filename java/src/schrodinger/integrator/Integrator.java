package schrodinger.integrator;

import schrodinger.QuantumState;
import schrodinger.potential.SchrodingerSystem;

public interface Integrator {

    double propagate(double[] psi, int n, QuantumState state, Direction direction);

    // TODO computePerturbativeCorrection should probably go into a separate interface. Consider refactoring.
    default double computePerturbativeCorrection(QuantumState state) {
        return 0;
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
