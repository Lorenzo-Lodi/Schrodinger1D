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

    //    Exponentially Fitted Numerov (EFN) with fixed peripheral (off-diagonal) coefficients,
    //    a.k.a. EFN with classical β=1/12 and fitted γ(ν)
    //    It should be method "S1" from Vanden Berghe et al., 2007
    public static Integrator getVignoli() {
        return new Vignoli();
    }


    //  Exponentially fitted Numerov S1 method (from Vanden Berghe et al., 2007)
    public static Integrator getExponentiallyFitted() {

        return new ExponentiallyFitted();
    }


}
