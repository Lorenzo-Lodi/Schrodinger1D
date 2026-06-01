package schrodinger;

public class Utils {

    public static String stepMsg(double actualOverMaximumRatio) {
        if (Double.isNaN(actualOverMaximumRatio) || actualOverMaximumRatio < 0) {
            return "ERROR!";
        }

        if (actualOverMaximumRatio >= 2.0) return "COARSE!    ";
        if (actualOverMaximumRatio >= 1.0) return "coarse!    ";
        if (actualOverMaximumRatio >= 0.5) return "sufficient  ";
        if (actualOverMaximumRatio >= 0.2) return "Good        ";
        return "Superfine   ";

    }

}
