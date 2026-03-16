package schrodinger.integrator;

import java.util.function.DoubleUnaryOperator;

public class IntegratorSpeedBenchmark {

    public static void time_benchmark(Integrator integrator) {
        int npoints = 1000000;

        double[] psi = new double[npoints];
        double[] currentPsiPrime = new double[1];
        for (int k = 0; k < 15; k++) {
            psi[k] = Math.cos(k);
            currentPsiPrime[0] = -Math.sin(k);
        }

        DoubleUnaryOperator q = i -> 1e-3 * (0.1111 + 0.2222 * i / npoints + 0.03333 / (npoints * npoints) * i * i);
        DoubleUnaryOperator q1 = i -> 1e-3 * (0.2222 / npoints + 2 * 0.3333 * i / (npoints * npoints));
        DoubleUnaryOperator q2 = i -> 1e-3 * (2 * 0.3333 / (npoints * npoints));
        long t0 = System.nanoTime();
        for (int n = 14; n < npoints - 1; n++) {
            psi[n + 1] = integrator.propagate(psi, currentPsiPrime, n, 0.123, q, q1, q2, Integrator.Direction.FORWARD);
        }
        long elapsedNanos = System.nanoTime() - t0;
        double timePerPoint = ((double) elapsedNanos) / (npoints);
        String msg = String.format("Time taken for integrating %d points is %10.3f ns / point for %s", npoints, timePerPoint, integrator.getClass().getSimpleName());
        System.out.println(msg);

//        Results 2026-03-05 h 16:40, commit e61c3cd0d1e0c9103dc4e39340a33adb5cf5c48f
//        Laptop i7 1355U, Amazon corretto 21.0.8
//        Time taken for integrating 1000000 points is     19.651 ns / point for EfnFixedBeta
//        Time taken for integrating 1000000 points is     36.674 ns / point for ExponentiallyFitted
//        Time taken for integrating 1000000 points is     29.706 ns / point for Numerov
//        Time taken for integrating 1000000 points is     60.352 ns / point for PredictorCorrector6
//        Time taken for integrating 1000000 points is     76.251 ns / point for PredictorCorrector8NumerovIter1
//        Time taken for integrating 1000000 points is    102.563 ns / point for PredictorCorrector8NumerovIter2
//        Time taken for integrating 1000000 points is     18.324 ns / point for Stormer5
//        Time taken for integrating 1000000 points is     14.718 ns / point for Stormer6
//        Time taken for integrating 1000000 points is     57.302 ns / point for Stormer8
//        Time taken for integrating 1000000 points is      9.685 ns / point for TaylorThreePoints

//        Results 2026-03-06 h 09:06, commit e77ee50e7a6e5cb8922cd0820b5503ea52a615fe
//        Desktop i5 4460, Oracle OpenJDK 19.0.1
//        Time taken for integrating 1000000 points is     60.920 ns / point for EfnFixedBeta
//        Time taken for integrating 1000000 points is     78.005 ns / point for ExponentiallyFitted
//        Time taken for integrating 1000000 points is    104.782 ns / point for Numerov
//        Time taken for integrating 1000000 points is    135.936 ns / point for PredictorCorrector6
//        Time taken for integrating 1000000 points is    226.429 ns / point for PredictorCorrector8NumerovIter1
//        Time taken for integrating 1000000 points is    315.115 ns / point for PredictorCorrector8NumerovIter2
//        Time taken for integrating 1000000 points is     33.766 ns / point for Stormer5
//        Time taken for integrating 1000000 points is     44.969 ns / point for Stormer6
//        Time taken for integrating 1000000 points is    101.723 ns / point for Stormer8
//        Time taken for integrating 1000000 points is     16.698 ns / point for TaylorThreePoints

//        Note 1: The i7 1355U is approx. 2.7 times faster than i5 4460, with a scatter of around 20% around this value.
//        These are the speed normalized to the TaylorThreePoints method:
//                                        i7-1355U  i5-4460
//        EfnFixedBeta	                     2.03	 3.65
//        ExponentiallyFitted	             3.79	 4.67
//        Numerov	                         3.07	 6.28
//        PredictorCorrector6	             6.23	 8.14
//        PredictorCorrector8NumerovIter1 	 7.87	13.56
//        PredictorCorrector8NumerovIter2	10.59	18.87
//        Stormer5	                         1.89	 2.02
//        Stormer6	                         1.52	 2.69
//        Stormer8	                         5.92	 6.09
//        TaylorThreePoints              	 1.00	 1.00

// These are the averages across the two machines, and sorted by speed:
//        PredictorCorrector8NumerovIter2	14.73
//        PredictorCorrector8NumerovIter1	10.72
//        PredictorCorrector6	             7.19
//        Stormer8	                         6.00
//        Numerov	                         4.67
//        ExponentiallyFitted	             4.23
//        EfnFixedBeta	                     2.84
//        Stormer6	                         2.11
//        Stormer5	                         1.96
//        TaylorThreePoints	                 1.00

// NOTE: TaylorThreePoints is very fast, which is no surprise as it's super simple, but its convergence is really poor.
// The real baseline is Numerov, and we can see that all methods apart from Predictor-Corrector have similar (usually BETTER)
//  speed! This really shouldn't be the case, as they are more complex.

    }

}
