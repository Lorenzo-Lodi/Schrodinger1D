package net.tinvention.schrodinger.grid;

public class GridFactory {

    private GridFactory() {
    }

    public static Grid generateUniformGrid(double rMin, double rMax, int numberOfPoints) {
        return new UniformGrid(rMin, rMax, numberOfPoints);
    }

    public static Grid generateLogarithmicGrid(double rMin, double rMax, int numberOfPoints, double rRef) {
        return new LogarithmicGrid(rMin, rMax, numberOfPoints, rRef);
    }

    public static Grid generateSurkusGrid(double rMin, double rMax, int numberOfPoints, double rRef, double alpha) {
        return new SurkusGrid(rMin, rMax, numberOfPoints, rRef, alpha);
    }

    public static Grid generateSqrtGrid(double rMin, double rMax, int numberOfPoints, double rRef) {
        return new SqrtGrid(rMin, rMax, numberOfPoints, rRef);
    }
}
