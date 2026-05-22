package schrodinger;

import java.util.function.DoubleUnaryOperator;

public class FractionalGridCache {
    private final double[] cache;
    private final int numFractions;
    private final double[] fractionalOffsets;
    public int nOfCacheHits = 0;
    public int nOfCacheMisses = 0;

    // Accepts any function taking in a (continuous) grid index i and tabularizes it for i=0, ... i = numGridPoints-1
    public FractionalGridCache(int numGridPoints, DoubleUnaryOperator function) {
        this(numGridPoints, new double[1], function);
    }

    public FractionalGridCache(int numGridPoints, double[] offsets, DoubleUnaryOperator function) {
        this.numFractions = offsets.length;
        this.fractionalOffsets = offsets;
        this.cache = new double[numGridPoints * numFractions];

        for (int i = 0; i < numGridPoints; i++) {
            for (int f = 0; f < numFractions; f++) {
                double continuousIndex = i + offsets[f];
                int flatIndex = (i * numFractions) + f;
                cache[flatIndex] = function.applyAsDouble(continuousIndex);
            }
        }
    }

    // O(numFractions) lookup
    public double get(double i) {
        return cache[computeCacheVectorIndex(i)];
    }

    public double getByFlatIndex(int flatIndex) {
        return cache[flatIndex];
    }

    public int computeCacheVectorIndex(double i) {
        int baseIndex = (int) i;
        double fraction = i - baseIndex;

        int fractionIndex = -1;
        for (int f = 0; f < numFractions; f++) {
            if (Math.abs(fraction - fractionalOffsets[f]) < 1e-9) {
                fractionIndex = f;
                break;
            }
        }
        if (fractionIndex >= 0) {
            nOfCacheHits++;
            return (baseIndex * numFractions) + fractionIndex;
        } else {
            nOfCacheMisses++;
            return -1;
        }
    }

    public int cacheVectorLength() {
        return cache.length;
    }
}
