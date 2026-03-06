package schrodinger;

import java.util.function.DoubleUnaryOperator;

public class FractionalGridCache {
    private final double[] cache;
    private final int numFractions;
    private final double[] fractionalOffsets;

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

    // O(1) fast lookup
    public double get(double i) {
        int baseIndex = (int) i;
        double fraction = i - baseIndex;

        int fractionIndex = 0;
        for (int f = 0; f < numFractions; f++) {
            if (Math.abs(fraction - fractionalOffsets[f]) < 1e-9) {
                fractionIndex = f;
                break;
            }
        }
        return cache[(baseIndex * numFractions) + fractionIndex];
    }
}
