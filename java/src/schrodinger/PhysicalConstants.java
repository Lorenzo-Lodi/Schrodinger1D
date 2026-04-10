package schrodinger;

public class PhysicalConstants {
    //    Source: 2022 CODATA, from the NIST Reference on Constants, Units and Uncertainty
    //    https://physics.nist.gov/cuu/Constants/index.html
    public static final double HARTREE_TO_INVERSE_CM = 219474.63136314;
    public static final double BOHR_TO_ANG = 0.529177210544;
    public static final double UMA_TO_ELECTRON_MASS = 1822.8884862781415;

    public static double toInverseCm(double energyHartree) {
        return energyHartree * HARTREE_TO_INVERSE_CM;
    }

    public static double toHartree(double energyInverseCm) {
        return energyInverseCm / HARTREE_TO_INVERSE_CM;
    }

}
