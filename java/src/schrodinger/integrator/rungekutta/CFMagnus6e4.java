package schrodinger.integrator.rungekutta;

import schrodinger.integrator.Integrator;

import java.util.function.DoubleUnaryOperator;

public final class CFMagnus6e4 extends CFMagnusAbstract implements Integrator {

    // 6th order requires 3 Gauss-Legendre nodes and 4 exponentials
    private static final int NODES = 3;
    private static final int EXPONENTIALS = 4;

    private static final double[] C = new double[NODES];

    private static final double[] A1W = new double[EXPONENTIALS];
    private static final double[][] W = new double[EXPONENTIALS][NODES];

    static {
        // 1. Gauss-Legendre Nodes mapped to [0, 1]
        double sqrt35 = Math.sqrt(3.0 / 5.0);
        C[0] = 0.5 - 0.5 * sqrt35;
        C[1] = 0.5;
        C[2] = 0.5 + 0.5 * sqrt35;

        // 2. Quadrature weights for interval size 1
        final double[] GW = new double[NODES];
        GW[0] = 5.0 / 18.0;
        GW[1] = 8.0 / 18.0;  // or 4.0 / 9.0
        GW[2] = 5.0 / 18.0;

        // 3. User-provided coefficients for CF6:4 (First half)
        // From Table 3 of 1102.5071v2.pdf
        // F arrays: 4 exponentials, 3 Legendre coefficients (A_1, A_2, A_3)
        final double[][] F = new double[EXPONENTIALS][3];
        F[0][0] = 1.0798524263824308825;
        F[0][1] = F[0][0] - (2.0 / 3.0) * F[0][0] * F[0][0];
        F[0][2] = 1.0 / (10.0 - 10.0 * F[0][0]);

        F[1][0] = 0.5 - F[0][0];
        F[1][1] = (1.0 - 4.0 * F[0][0] + 2.0 * F[0][0] * F[0][0]) / 3.0;
        F[1][2] = -F[0][2];

        // 4. Apply Time-Reversal Symmetry for the second half
        // F[2] is the mirror of F[1], F[3] is the mirror of F[0]
        for (int i = 0; i < 2; i++) {
            F[3 - i][0] = F[i][0]; // Even parity
            F[3 - i][1] = -F[i][1]; // Odd parity flips
            F[3 - i][2] = F[i][2]; // Even parity
        }

        // 5. Precompute effective node weights using exact polynomials
        for (int i = 0; i < EXPONENTIALS; i++) {
            A1W[i] = F[i][0];

            for (int k = 0; k < NODES; k++) {
                double c = C[k];

                // shifted Legendre polynomials on [0,1]
                double p0 = 1.0;
                double p1 = 2.0 * c - 1.0;
                double p2 = 6.0 * c * c - 6.0 * c + 1.0;


                // Eq. (25): A_n has factor (2n - 1) -> 1, 3, 5, 7
                W[i][k] = GW[k] * (
                        F[i][0] * p0
                                + 3.0 * F[i][1] * p1
                                + 5.0 * F[i][2] * p2
                );
            }
        }

    }

    public CFMagnus6e4() {
        super(NODES, EXPONENTIALS, C, W, A1W);
    }

    @Override
    public int minHistoryLength() {
        return 1;
    }

    @Override
    public int globalConvergenceOrder() {
        return 6;
    }

}









