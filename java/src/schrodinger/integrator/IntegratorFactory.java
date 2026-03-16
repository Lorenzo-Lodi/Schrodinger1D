package schrodinger.integrator;

import schrodinger.integrator.numerovlike.EFN;
import schrodinger.integrator.numerovlike.Numerov;
import schrodinger.integrator.numerovlike.EFNFixedBeta;
import schrodinger.integrator.predcorr.PredictorCorrector6;
import schrodinger.integrator.predcorr.PredictorCorrector8i1;
import schrodinger.integrator.predcorr.PredictorCorrector8i2;
import schrodinger.integrator.rungekutta.RK45DP;
import schrodinger.integrator.rungekutta.RKN4;
import schrodinger.integrator.stormer.Stormer5;
import schrodinger.integrator.stormer.Stormer6;
import schrodinger.integrator.stormer.Stormer8;
import schrodinger.integrator.stormer.Stormer8i;
import schrodinger.integrator.various.Obrechkoff6;
import schrodinger.integrator.various.TaylorThreePoints;

import java.util.ArrayList;
import java.util.List;

public class IntegratorFactory {

    public static Integrator getTaylorThreePoints() {
        return new TaylorThreePoints();
    }

    public static Integrator getNumerov() {
        return new Numerov();
    }

    public static Integrator getEFNFixedBeta() {
        return new EFNFixedBeta();
    }

    public static Integrator getEFN() {
        return new EFN();
    }

    public static Integrator getStormer5() {
        return new Stormer5();
    }

    public static Integrator getStormer6() {
        return new Stormer6();
    }

    public static Integrator getStormer8() {
        return new Stormer8();
    }

    public static Integrator getStormer8i() {
        return new Stormer8i();
    }

    public static Integrator getPredictorCorrector6() {
        return new PredictorCorrector6();
    }

    public static Integrator getPredictorCorrector8i1() {
        return new PredictorCorrector8i1();
    }

    public static Integrator getPredictorCorrector8i2() {
        return new PredictorCorrector8i2();
    }

    public static Integrator getObrechkoff6() {
        return new Obrechkoff6();
    }

    public static Integrator getRKN4() {
        return new RKN4();
    }

    public static Integrator getRK45DP() {
        return new RK45DP();
    }

    public static List<Integrator> getAll() {
        List<Integrator> list = new ArrayList<>();
        list.add(getTaylorThreePoints());
        list.add(getNumerov());
        list.add(getEFNFixedBeta());
        list.add(getEFN());
        list.add(getStormer5());
        list.add(getStormer6());
        list.add(getStormer8());
        list.add(getStormer8i());
        list.add(getObrechkoff6());
        list.add(getPredictorCorrector6());
        list.add(getPredictorCorrector8i1());
        list.add(getPredictorCorrector8i2());
        list.add(getRKN4());
        list.add(getRK45DP());
        return list;
    }

}