package schrodinger.integrator.verification;

public class ConvergenceParams {
    public final double convergenceOrder;
    public final double characteristicStepSize;

    public ConvergenceParams(double convergenceOrder, double characteristicStepSize) {
        this.convergenceOrder = convergenceOrder;
        this.characteristicStepSize = characteristicStepSize;
    }

    @Deprecated
    public double calculateMaxError(int np) {
        return 1e-10;
    }

    public double calculateMaxError(double stepSize) {
        return Math.pow(stepSize / characteristicStepSize, convergenceOrder);
    }
}
