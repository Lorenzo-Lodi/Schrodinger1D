package schrodinger.integrator;

public class IntegratorFactory {

    private IntegratorFactory() {
    }

    public static Integrator getTaylorThreePoints() {
        return new TaylorThreePoints();
    }

    public static Integrator getNumerov() {
        return new Numerov();
    }

    public static Integrator getVignoli() {
        return new Vignoli();
    }

}
