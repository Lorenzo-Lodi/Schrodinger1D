package schrodinger.grid;

public class Grid {
    private final double rMin;
    private final double rMax;
    private final double yMin;
    private final double yMax;
    private final double stepSizeYCoordinate;
    private final int numberOfPoints;
    private final MappingStrategy mappingStrategy;

    public Grid(double rMin, double rMax, int numberOfPoints, MappingStrategy mappingStrategy) {
        this.rMin = Math.min(rMin, rMax);
        this.rMax = Math.max(rMin, rMax);
        this.numberOfPoints = Math.max(numberOfPoints, 2);
        this.mappingStrategy = mappingStrategy;

        this.yMax = this.y(this.rMax);
        this.yMin = this.y(this.rMin);
        this.stepSizeYCoordinate = (yMax - yMin) / (numberOfPoints - 1);
    }

    public double getFirstYValue() {
        return yAtGridPoint(0);
    }

    public double getLastYValue() {
        return yMax;
    }

    public double getStepSizeYCoordinate() {
        return stepSizeYCoordinate;
    }

    public int getNumberOfPoints() {
        return numberOfPoints;
    }

    /**
     *
     * @param i Index going from 0 to numberOfPoints-1
     * @return Value of the i-th grid point
     */
    public double yAtGridPoint(double i) {
        return yMin + getStepSizeYCoordinate() * i;
    }

    public double rAtGridPoint(double i) {
        return r(yAtGridPoint(i));
    }

    /**
     * The mapping function y(r), mapping the r coordinate to the transformed y
     * coordinate
     *
     * @param r Value of the r coordinate
     * @return Corresponding value y(r)
     */
    public double y(double r) {
        return mappingStrategy.y(r);
    }

    /**
     * The inverse mapping function r(y), mapping the transformed y coordinate to r
     *
     * @param y Value of the mapped coordinate
     * @return Corresponding value r(y)
     */
    public double r(double y) {
        return mappingStrategy.r(y);
    }

    /**
     * It's the derivative of the r(y) function, referred to as g(y) in Meshkov2008
     *
     * @param y Value of the mapped coordinate
     * @return Value of g(y)=r'(y)
     */
    public double g(double y) {
        return mappingStrategy.g(y);
    }

    /**
     * It's the function F(y) of eq. (9) in Meshkov2008, the additional term to the
     * transformed function \tilde{Q}
     *
     * @param y
     * @return
     */
    public double F(double y) {
        return mappingStrategy.F(y);
    }

    public double gAtGridPoint(double i) {
        return g(yAtGridPoint(i));
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("i           R                Y                g(y)             F(y)\n");
        for (int i = 0; i < getNumberOfPoints(); i++) {
            sb.append(String.format("%8d", i));
            sb.append("   ");
            double y = yAtGridPoint(i);
            double r = r(y);
            if (r >= 0.d) {
                sb.append(" ");
            }
            sb.append(String.format("%.012f", r) + "  ");
            if (y >= 0.d) {
                sb.append(" ");
            }
            sb.append(String.format("%.012f", y));
            sb.append("  ");
            sb.append(String.format("%.012f", g(y)));
            sb.append("  ");
            sb.append(String.format("%.012f", F(y)) + "\n");
        }
        return sb.toString();
    }

}
