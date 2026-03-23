package schrodinger.integrator.rungekutta;

import schrodinger.integrator.Integrator;

import java.util.function.DoubleUnaryOperator;

// F.Casas, S.Blanes, A.Escorihuela-Tomàs, Runge-Kutta-Nyström symplectic splitting methods of order 8
// https://arxiv.org/abs/2202.01541

//The methods in the Blanes paper (A17, B17, etc.) are highly optimized specifically for autonomous systems: $y''(x) = f(y)$.
//For an autonomous system, certain Lie-algebraic commutators naturally evaluate to zero. Blanes et al. took advantage of this to drastically reduce the number of equations they had to solve to reach order 8.
//However, your equation is non-autonomous: $y''(x) = -Q(x)y(x)$. Because the force depends explicitly on $x$, the time derivatives like $Q'(x)$ and $Q''(x)$ generate new non-zero terms in the Taylor expansion.

//In short: it gives only h^4 convergence instead of h^8, so it's unusuitable!

public final class RKN8Blanes implements Integrator {

    // 17-stage symmetric splitting methods B17
    private static final int STAGES = 17;

    private static final double[] A = new double[STAGES];
    private static final double[] B = new double[STAGES + 1];
    private static final double[] C = new double[STAGES];

    static {
        double a1 = 0.160227696073839513690970240076;
        double a2 = 0.306354507436867319879440957100;
        double a3 = 0.308395508895171191756544975556;
        double a4 = 0.120362086566233408450063177659;
        double a5 = -0.622888687549183872072186218718;
        double a6 = 0.635560951632990078378672016548;
        double a7 = -0.144226974795419229640437363913;
        double a8 = -0.284867527074173816678992817545;

        double b1 = 0.0514196142537210073343152693459;
        double b2 = 0.250497030318342871458417941091;
        double b3 = 0.512412268300327350035492806653;
        double b4 = -0.231597138650894401279645184364;
        double b5 = 0.116091323536875759881216298975;
        double b6 = -0.0098365173246965763985763034283;
        double b7 = -0.108032771466281638634277563747;
        double b8 = 0.249039864198023642002940910070;

        // The center weights are constrained to guarantee time-consistency (sum to 1)
        double sumA = a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8;
        double a9 = 1.0 - 2.0 * sumA;

        double sumB = b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8;
        double b9 = 0.5 - sumB;

        // 1. Build the first half of the composition
        A[0] = a1;
        B[0] = b1;
        A[1] = a2;
        B[1] = b2;
        A[2] = a3;
        B[2] = b3;
        A[3] = a4;
        B[3] = b4;
        A[4] = a5;
        B[4] = b5;
        A[5] = a6;
        B[5] = b6;
        A[6] = a7;
        B[6] = b7;
        A[7] = a8;
        B[7] = b8;

        // Center node
        A[8] = a9;
        B[8] = b9;

        // 2. Apply symmetry for the second half: b_{m+1-i} = b_i and a_{m-i} = a_i
        for (int i = 0; i < 8; i++) {
            A[16 - i] = A[i];     // a_17 to a_10
            B[17 - i] = B[i];     // b_18 to b_11
        }
        B[9] = b9; // b_10 = b_9, maintaining the twin center splits

        // 3. Compute the C nodes (cumulative sum of the B drift steps)
        double currentC = 0.0;
        for (int i = 0; i < STAGES; i++) {
            currentC += B[i];
            C[i] = currentC;
        }
    }

    @Override
    public int minHistoryLength() {
        return 1;
    }

    @Override
    public int globalConvergenceOrder() {
        return 8;
    }

    @Override
    public double propagate(double[] psi, double[] currentPsiPrime, int n, double step,
                            DoubleUnaryOperator qTilde,
                            DoubleUnaryOperator qTildePrime,
                            DoubleUnaryOperator qTildeDoublePrime,
                            Direction direction) {

        double h = step;
        int d = direction.getValue();

        // Initialize iteration with the current step's states
        double y = psi[n];
        double yp = currentPsiPrime[0];

        // Process the 17-stage sequence alternating Position and Velocity
        for (int i = 0; i < STAGES; i++) {

            // 1. Position Drift
            y += h * d * B[i] * yp;

            // 2. Evaluate force function at the new intermediate time node C_i
            double q = qTilde.applyAsDouble(n + C[i] * d);
            double f = -q * y;

            // 3. Velocity Kick
            yp += h * d * A[i] * f;
        }

        // Final position drift to conclude the time step
        y += h * d * B[STAGES] * yp;

        // Update momentum reference and return the newly calculated position
        currentPsiPrime[0] = yp;
        return y;
    }
}
