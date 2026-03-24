package schrodinger.integrator.rungekutta;

import schrodinger.integrator.Integrator;

import java.util.function.DoubleUnaryOperator;

public final class CFMagnus8 extends CFMagnusAbstract implements Integrator {

    private static final int NODES = 4;
    private static final int EXPONENTIALS = 11;

    private static final double[] C = new double[NODES];

    private static final double[] A1W = new double[EXPONENTIALS];
    private static final double[][] W = new double[EXPONENTIALS][NODES];

    static {
        double sqrt30 = Math.sqrt(30.0);

        // Nodes mapped to [0, 1] for the user's qTilde function
        C[0] = 0.5 - Math.sqrt((15.0 + 2.0 * sqrt30) / 140.0);
        C[1] = 0.5 - Math.sqrt((15.0 - 2.0 * sqrt30) / 140.0);
        C[2] = 0.5 + Math.sqrt((15.0 - 2.0 * sqrt30) / 140.0);
        C[3] = 0.5 + Math.sqrt((15.0 + 2.0 * sqrt30) / 140.0);

        // Quadrature weights for interval size 1
        final double[] GW = new double[NODES];
        GW[0] = (18.0 - sqrt30) / 72.0;
        GW[1] = (18.0 + sqrt30) / 72.0;
        GW[2] = (18.0 + sqrt30) / 72.0;
        GW[3] = (18.0 - sqrt30) / 72.0;

        final double[][] F = new double[EXPONENTIALS][4];
        F[0][0] = 0.169715531043933180094151;
        F[0][1] = 0.152866146944615909929839;
        F[0][2] = 0.119167378745981369601216;
        F[0][3] = 0.068619226448029559107538;

        F[1][0] = 0.379420807516005431504230;
        F[1][1] = 0.148839980923180990943008;
        F[1][2] = -0.115880829186628075021088;
        F[1][3] = -0.188555246668412628269760;

        F[2][0] = 0.469459306644050573017994;
        F[2][1] = -0.379844237839363505173921;
        F[2][2] = 0.022898814729462898505141;
        F[2][3] = 0.571855043580130805495594;

        F[3][0] = -0.448225927391070886302766;
        F[3][1] = 0.362889857410989942809900;
        F[3][2] = -0.022565582830528472333301;
        F[3][3] = -0.544507517141613383517695;

        F[4][0] = -0.293924473106317605373923;
        F[4][1] = -0.026255628265819381983204;
        F[4][2] = 0.096761509131620390100068;
        F[4][3] = 0.000018330145571671744069;

        F[5][0] = 0.447109510586798614120629;
        F[5][1] = 0.0;
        F[5][2] = -0.200762581179816221704073;
        F[5][3] = 0.0;

        // Apply Time-Reversal Symmetry
        for (int i = 0; i < 5; i++) {
            F[10 - i][0] = F[i][0];
            F[10 - i][1] = -F[i][1];
            F[10 - i][2] = F[i][2];
            F[10 - i][3] = -F[i][3];
        }

        // Precompute effective node weights from Table 4
        for (int i = 0; i < EXPONENTIALS; i++) {
            A1W[i] = F[i][0];

            for (int k = 0; k < NODES; k++) {
                double x = C[k];

                // shifted Legendre polynomials on [0,1]
                double p0 = 1.0;
                double p1 = 2.0 * x - 1.0;
                double p2 = 6.0 * x * x - 6.0 * x + 1.0;
                double p3 = 20.0 * x * x * x - 30.0 * x * x + 12.0 * x - 1.0;

                // Eq. (25): A_n has factor (2n-1)
                W[i][k] = GW[k] * (
                        F[i][0] * p0
                                + 3.0 * F[i][1] * p1
                                + 5.0 * F[i][2] * p2
                                + 7.0 * F[i][3] * p3
                );
            }
        }

    }

    public CFMagnus8() {
        super(NODES, EXPONENTIALS, C, W, A1W);
    }

    @Override
    public int minHistoryLength() {
        return 1;
    }

    @Override
    public int globalConvergenceOrder() {
        return 8;
    }

}
