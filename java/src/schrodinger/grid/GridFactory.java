package schrodinger.grid;

public class GridFactory {

    private GridFactory() {
    }

    public static Grid generateUniformGrid(double rMin, double rMax, int numberOfPoints) {
        return new Grid(rMin, rMax, numberOfPoints, new MappingStrategy() {
        });
    }

    public static Grid generateLogarithmicGrid(double rMin, double rMax, int numberOfPoints, double rRef) {
        return new Grid(rMin, rMax, numberOfPoints, new MappingLogarithmic(rRef));
    }

    public static Grid generateSurkusGrid(double rMin, double rMax, int numberOfPoints, double rRef, double alpha) {
        return new Grid(rMin, rMax, numberOfPoints, new MappingSurkus(rRef, alpha));
    }

    public static Grid generateSqrtGrid(double rMin, double rMax, int numberOfPoints, double rRef) {
        return new Grid(rMin, rMax, numberOfPoints, new MappingSqrt(rRef));
    }
}
