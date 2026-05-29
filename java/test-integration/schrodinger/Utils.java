package schrodinger;

public class Utils {

    public static String stepMsg(double actualOverMaximumRatio) {
        if (actualOverMaximumRatio > 2.) {
            return "!!";
        }
        if (actualOverMaximumRatio > 1.) {
            return "!";
        }
        if (actualOverMaximumRatio > 0.5) {
            return "careful";
        }
        if (actualOverMaximumRatio > 0.2) {
            return "OK";
        }
        return "Superfine";
    }

}
