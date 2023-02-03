package net.tinvention.schrodinger.grid;

public class UniformGrid implements Grid {
	private double rMin;
	private double rMax;
	private double stepSize;
	private int numberOfPoints;

	public UniformGrid(double rMin, double rMax, int numberOfPoints) {
		this.rMin = Math.min(rMin, rMax);
		this.rMax = Math.max(rMin, rMax);
		this.numberOfPoints = Math.max(numberOfPoints, 2);
		this.stepSize = (rMax - rMin) / (numberOfPoints - 1);
	}

	public UniformGrid(double rMin, double rMax, double desiredStepSize) {
		this(rMin, rMax, (int) (Math.abs(rMax - rMin) / desiredStepSize));
	}

	@Override
	public double getFirstYValue() {
		return rMin;
	}

	@Override
	public double getLastYValue() {
		return rMax;
	}

	@Override
	public double getStepSizeYCoordinate() {
		return stepSize;
	}

	@Override
	public int getNumberOfPoints() {
		return numberOfPoints;
	}

	@Override
	public double mappingFunctionYofR(double r) {
		return r;
	}

	@Override
	public double mappingFunctionRofY(double y) {
		return y;
	}

	@Override
	public double mappingFunctionGofY(double y) {
		return 1.d;
	}

	@Override
	public double mappingFunctionFofY(double y) {
		return 0;
	}

}
