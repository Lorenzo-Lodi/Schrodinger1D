package net.tinvention.schrodinger.grid;

public class SurkusGrid implements Grid {
	private double rMin;
	private double rMax;
	private double stepSize;
	private int numberOfPoints;
	private double rRef;
	private double alpha;
	private double yMin;
	private double yMax;

	public SurkusGrid(double rMin, double rMax, int numberOfPoints, double rRef, double alpha) {
		this.rMin = Math.min(rMin, rMax);
		this.rMax = Math.max(rMin, rMax);
		this.rRef = rRef;
		this.alpha = alpha;
		this.numberOfPoints = Math.max(numberOfPoints, 2);
		this.yMax = this.mappingFunctionYofR(this.rMax);
		this.yMin = this.mappingFunctionYofR(this.rMin);
		this.stepSize = (yMax - yMin) / (numberOfPoints - 1);
	}

	@Override
	public double getFirstYValue() {
		return this.mappingFunctionYofR(rMin);
	}

	@Override
	public double getLastYValue() {
		return this.mappingFunctionYofR(rMax);
	}

	@Override
	public double getStepSizeYCoordinate() {
		return this.stepSize;
	}

	@Override
	public int getNumberOfPoints() {
		return this.numberOfPoints;
	}

	@Override
	public double mappingFunctionYofR(double r) {
		return (Math.pow(r / rRef, alpha) - 1.0d) / (Math.pow(r / rRef, alpha) + 1.0d);
	}

	@Override
	public double mappingFunctionRofY(double y) {
		return rRef * Math.pow((1 + y) / (1 - y), 1.0d / alpha);
	}

	@Override
	public double mappingFunctionGofY(double y) {
		return (2.0d * rRef * Math.pow(1 + y, -1.0d + 1.0d / alpha))
				/ (alpha * Math.pow(1.0d - y, 1.0d + 1.0d / alpha));
	}

	@Override
	public double mappingFunctionFofY(double y) {
		return (1.0d - 1.0d / Math.pow(alpha, 2)) / Math.pow(1 - y * y, 2);
	}

}
