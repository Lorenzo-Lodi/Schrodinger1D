package schrodinger.integrator;

import schrodinger.integrator.multi_step.numerovlike.EFN;
import schrodinger.integrator.multi_step.numerovlike.Numerov;
import schrodinger.integrator.multi_step.numerovlike.EFNFixedBeta;
import schrodinger.integrator.multi_step.predcorr.PC6;
import schrodinger.integrator.multi_step.predcorr.PC8i1;
import schrodinger.integrator.multi_step.predcorr.PC8i2;
import schrodinger.integrator.one_step.*;
import schrodinger.integrator.multi_step.Cowell5;
import schrodinger.integrator.multi_step.Cowell6;
import schrodinger.integrator.multi_step.Stormer8;
import schrodinger.integrator.multi_step.Cowell8;
import schrodinger.integrator.multi_step.Obrechkoff6;
import schrodinger.integrator.multi_step.Verlet;

import java.util.ArrayList;
import java.util.List;

public class IntegratorFactory {

    private IntegratorFactory() {
    }

    public static Integrator getVerlet() {
        return new Verlet();
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

    public static Integrator getCowell5() {
        return new Cowell5();
    }

    public static Integrator getCowell6() {
        return new Cowell6();
    }

    public static Integrator getStormer8() {
        return new Stormer8();
    }

    public static Integrator getCowell8() {
        return new Cowell8();
    }

    public static Integrator getPC6() {
        return new PC6();
    }

    public static Integrator getPC8i1() {
        return new PC8i1();
    }

    public static Integrator getPC8i2() {
        return new PC8i2();
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

    public static Integrator getCFMagnus4() {
        return new CFMagnus4();
    }

    public static Integrator getCFMagnus6e4() {
        return new CFMagnus6e4();
    }

    public static Integrator getCFMagnus6e5Opt() {
        return new CFMagnus6e5Opt();
    }

    public static Integrator getCFMagnus8() {
        return new CFMagnus8();
    }

    public static Integrator getBestOneStepIntegrator() {
        return new CFMagnus8();
    }


    public static List<Integrator> getAll() {
        List<Integrator> list = new ArrayList<>();
        list.add(getVerlet());
        list.add(getNumerov());
        list.add(getEFNFixedBeta());
        list.add(getEFN());
        list.add(getCowell5());
        list.add(getCowell6());
        list.add(getStormer8());
        list.add(getCowell8());
        list.add(getObrechkoff6());
        list.add(getPC6());
        list.add(getPC8i1());
        list.add(getPC8i2());
        list.add(getRKN4());
        list.add(getRK45DP());
        list.add(getCFMagnus4());
        list.add(getCFMagnus6e4());
        list.add(getCFMagnus6e5Opt());
        list.add(getCFMagnus8());
        return list;
    }

}