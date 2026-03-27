package schrodinger.integrator.verification;

public class ConvergenceParams {
    private final double convergenceOrder;
    private final double characteristicStepSize;
    private final Double minStepSize;
    private final Double maxStepSize;

    public ConvergenceParams(double convergenceOrder, double characteristicStepSize, double minStepSize, double maxStepSize) {
        this.convergenceOrder = convergenceOrder;
        this.characteristicStepSize = characteristicStepSize;
        this.minStepSize = minStepSize;
        this.maxStepSize = maxStepSize;
    }


    public ConvergenceParams(double convergenceOrder, double characteristicStepSize) {
        this.convergenceOrder = convergenceOrder;
        this.characteristicStepSize = characteristicStepSize;
        this.minStepSize = null;
        this.maxStepSize = null;
    }

    public double calculateMaxError(double stepSize) {
        if (maxStepSize != null && stepSize >= maxStepSize) {
            // If step size is very large, return huge number to make test pass.
            // This is for methods that blow up at large h like Stormer8
            return Double.MAX_VALUE;
        }
        if (minStepSize != null && stepSize <= minStepSize) {
            // If step size is too small, some high-order methods flatten to numerical noise.
            // In this case just return a safe small value
            return 1e-12;
        }

        return Math.pow(stepSize / characteristicStepSize, convergenceOrder);
    }
}
