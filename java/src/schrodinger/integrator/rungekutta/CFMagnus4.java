package schrodinger.integrator.rungekutta;

import schrodinger.integrator.Integrator;

public final class CFMagnus4 extends CFMagnusAbstract implements Integrator {

    private static final int NODES = 2;
    private static final int EXPONENTIALS = 2;
    private static final double[] C = new double[NODES]; // Gauss-Legendre nodes
    private static final double[] A1W = new double[EXPONENTIALS];
    private static final double[][] W = new double[EXPONENTIALS][NODES]; // Weights

    static {
        double sqrt3 = Math.sqrt(3.0);

        // 1. Gauss-Legendre Nodes on [0, 1]
        C[0] = 0.5 - sqrt3 / 6.0;
        C[1] = 0.5 + sqrt3 / 6.0;

        // 2. Weights for Exponential 1 (assumes reverse iteration)
        W[1][0] = 0.25 + sqrt3 / 6.0;
        W[1][1] = 0.25 - sqrt3 / 6.0;

        // 3. Weights for Exponential 2 (assumes reverse iteration)
        W[0][0] = 0.25 - sqrt3 / 6.0;
        W[0][1] = 0.25 + sqrt3 / 6.0;

        // 4. Precompute A1W (Sum of weights for each exponential)
        A1W[0] = W[0][0] + W[0][1]; // 0.5
        A1W[1] = W[1][0] + W[1][1]; // 0.5
    }

    public CFMagnus4() {
        super(NODES, EXPONENTIALS, C, W, A1W);
    }

    @Override
    public int minHistoryLength() {
        return 1;
    }

    @Override
    public int globalConvergenceOrder() {
        return 4;
    }

}
