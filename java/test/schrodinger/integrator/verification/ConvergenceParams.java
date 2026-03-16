package schrodinger.integrator.verification;

public class ConvergenceParams {
    public final double characteristicNumberOfPoints;
    public final double convergenceOrder;

    public ConvergenceParams(double convergenceOrder, double characteristicNumberOfPoints) {
        this.convergenceOrder = convergenceOrder;
        this.characteristicNumberOfPoints = characteristicNumberOfPoints;
    }

    public double calculateMaxError(int np) {
        return Math.pow(characteristicNumberOfPoints / np, convergenceOrder);
    }
}
