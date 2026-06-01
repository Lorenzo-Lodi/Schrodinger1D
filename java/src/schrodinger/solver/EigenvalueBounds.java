package schrodinger.solver;

import java.util.Arrays;

public class EigenvalueBounds {

    private final double[] lowerBounds;
    private final double[] upperBounds;
    private final int[] lowerNodes;
    private final int[] upperNodes;
    private static final double LARGE_VALUE = Double.MAX_VALUE / 16.0;

    public EigenvalueBounds(int maxQuantumNumber) {
        int size = maxQuantumNumber + 1;
        this.lowerBounds = new double[size];
        this.upperBounds = new double[size];
        this.lowerNodes = new int[size];
        this.upperNodes = new int[size];

        // Initialize to extreme limits
        Arrays.fill(lowerBounds, -LARGE_VALUE);
        Arrays.fill(upperBounds, LARGE_VALUE);
        Arrays.fill(lowerNodes, -1);
        Arrays.fill(upperNodes, -1);
    }

    public void updateBounds(double energy, int nNodes) {

        for (int v = 0; v < Math.min(nNodes, upperBounds.length); v++) {
            if (upperBounds[v] > energy) {
                upperBounds[v] = energy;
                upperNodes[v] = nNodes;
            }
        }

        for (int v = nNodes; v < lowerBounds.length; v++) {
            if (lowerBounds[v] < energy) {
                lowerBounds[v] = energy;
                lowerNodes[v] = v;
            }
        }

    }

    public boolean isUpperBoundDefined(int v) {
        if (v >= upperBounds.length) {
            return false;
        }
        return upperBounds[v] < LARGE_VALUE;
    }

    public boolean isLowerBoundDefined(int v) {
        if (v >= upperBounds.length) {
            return false;
        }
        return lowerBounds[v] > -LARGE_VALUE;
    }

    public double getUpperBound(int i) {
        return upperBounds[i];
    }

    public double getLowerBound(int i) {
        return lowerBounds[i];
    }

    public int getLowerNodes(int i) {
        return lowerNodes[i];
    }

    public int getUpperNodes(int i) {
        return upperNodes[i];
    }

}
