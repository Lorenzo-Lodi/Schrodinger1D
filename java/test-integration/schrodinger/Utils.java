package schrodinger;

public class Utils {

    public static String stepMsg(double actualOverMaximumRatio) {
        if (Double.isNaN(actualOverMaximumRatio) || actualOverMaximumRatio < 0) {
            return "ERROR!";
        }

        if (actualOverMaximumRatio >= 2.0) return "WARNING!    ";
        if (actualOverMaximumRatio >= 1.0) return "warning!    ";
        if (actualOverMaximumRatio >= 0.5) return "careful     ";
        if (actualOverMaximumRatio >= 0.2) return "OK          ";
        return "Superfine   ";

    }

}
