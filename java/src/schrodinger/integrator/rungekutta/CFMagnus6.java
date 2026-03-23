package schrodinger.integrator.rungekutta;

import schrodinger.integrator.Integrator;

import java.util.function.DoubleUnaryOperator;

public final class CFMagnus6 {

    static {
        // From Table 3 of 1102.5071v2.pdf
        // Four exponentials CF6:4
        double[][] F = new double[2][2];
        F[0][0] = 1.0798524263824308825; //0.5 + (5400 - 600*Sqrt(6))**0.3333333333333333/60. + ((9 + Sqrt(6))/5.)**0.3333333333333333/(2.*3**0.6666666666666666)
        F[0][1] = F[0][0] - (2. / 3.) * F[0][0] * F[0][0];
        F[0][2] = 1. / (10. - 10. * F[0][0]);

        F[1][0] = 0.5 - F[0][0];
        F[1][1] = (1. - 4. * F[0][0] + 2. * F[0][0] * F[0][0]) / 3.;
        F[1][2] = -F[0][2];

        F[2][0] = 0.;
        F[2][1] = 0.;
        F[2][2] = 0.;

        // Five exponentials CF6:5
        F[0][0] = 0.16;
        F[0][1] = 0.14587456942714338561;
//        F[1][0] =
//                F[1][1] =
//                        F[2][0] =
//                                F[2][1] =

    }

    public static void main(String[] args) {

    }

}
