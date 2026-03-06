package schrodinger;

import org.junit.jupiter.api.Test;

import java.util.function.DoubleUnaryOperator;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class FractionalGridCacheTest {

    @Test
    public void fractionalGridCacheTestBasicTestManualOffsets() {
        DoubleUnaryOperator function = i -> i * i + 2 * i + 3;
        int npoints = 10;
        double[] offsets = new double[1];
        offsets[0] = 0.0;
        FractionalGridCache cache = new FractionalGridCache(npoints, offsets, function);

        for (int i = 0; i < npoints; i++) {
            double val = cache.get(i);
            assertEquals(function.applyAsDouble(i), val);
        }
    }

    @Test
    public void fractionalGridCacheTestBasicTestNoOffsets() {
        DoubleUnaryOperator function = i -> i * i + 2 * i + 3;
        int npoints = 10;
        FractionalGridCache cache = new FractionalGridCache(npoints, function);

        for (int i = 0; i < npoints; i++) {
            double val = cache.get(i);
            assertEquals(function.applyAsDouble(i), val);
        }
    }

    @Test
    public void fractionalGridCacheTestBasicTestWithOffsets() {
        DoubleUnaryOperator function = i -> i * i + 2 * i + 3;
        int npoints = 10;
        double[] offsets = {0., 1. / 5., 1. / 3., 1. / 20.};
        FractionalGridCache cache = new FractionalGridCache(npoints, offsets, function);

        for (int i = 0; i < npoints; i++) {
            for (int f = 0; f < offsets.length; f++) {
                double continuousIndex = i + offsets[f];
                double val = cache.get(continuousIndex);
                assertEquals(function.applyAsDouble(continuousIndex), val);
            }
        }
    }


    @Test
    public void fractionalGridCacheTestPerformanceTest() {
        int npoints = 10000000;
        DoubleUnaryOperator function = i -> Math.sin(i / npoints) * Math.cos(i / npoints) + 2 * Math.pow(i / npoints, 0.123)
                + Math.log10(i) + 1. / Math.tan(i / npoints) + 3;
        double[] offsets = new double[1];
        offsets[0] = 0.0;
        FractionalGridCache cache = new FractionalGridCache(npoints, offsets, function);

        long t0 = System.nanoTime();
        double sumCached = 0.;
        for (int i = 0; i < npoints; i++) {
            sumCached += cache.get(i);
        }
        double elapsedNanosPerPointCached = ((double) (System.nanoTime() - t0)) / npoints;

        t0 = System.nanoTime();
        double sumNotCached = 0.;
        for (int i = 0; i < npoints; i++) {
            sumNotCached += function.applyAsDouble(i);
        }
        double elapsedNanosPerPointNotCached = ((double) (System.nanoTime() - t0)) / npoints;

        assertEquals(sumNotCached, sumCached);
        double ratioNotCachedOverCached = elapsedNanosPerPointNotCached / elapsedNanosPerPointCached;
//        System.out.println(npoints + " " + ratioNotCachedOverCached);
        assertTrue(ratioNotCachedOverCached > 10.);
    }
}


