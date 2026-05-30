package schrodinger.solver;

public class Bounds {

    private final Bound[] bounds;

    public Bounds(int vMax) {
        bounds = new Bound[vMax + 1];
        for (int v = 0; v < bounds.length; v++) {
            bounds[v] = new Bound();
        }
    }

    public int updateBounds(double energy, int nNodes) {
        int updatedBounds = 0;
        for (int v = 0; v < Math.min(nNodes, bounds.length); v++) {
            if (bounds[v].upperBound > energy) {
                bounds[v].upperBound = energy;
                updatedBounds++;
            }
        }

        for (int v = nNodes; v < bounds.length; v++) {
            if (bounds[v].lowerBound < energy) {
                bounds[v].lowerBound = energy;
                updatedBounds++;
            }
        }

        return updatedBounds;
    }

    private static class Bound {
        double lowerBound = -Double.MAX_VALUE / 16.;
        double upperBound = Double.MAX_VALUE / 16.;
    }

}
