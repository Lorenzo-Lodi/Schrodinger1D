package morse;

import schrodinger.QuantumLevel;
import schrodinger.SchrodingerSystem;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.IntegratorFactory;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialMorse;
import schrodinger.solver.RefinementStrategy;
import schrodinger.solver.ShootingSolver;

import static schrodinger.PhysicalConstants.*;

public class MorseMain {

    public static void main(String[] args) {

        double a = 1.5;
        double rmin = 2.;
        double De = toHartree(50000.);

        PhysicalPotential potential = new PhysicalPotentialMorse(rmin, a, De);
        double mass = 100. * UMA_TO_ELECTRON_MASS;

        double omega0 = a * Math.sqrt(2.0 * De / mass);
        double xe = omega0 / (4. * De);
        double A = 1. / xe;
        int nOfBoundStates = (int) ((A + 1.) * 0.5);
        System.out.printf("omega0=%20.8f cm-1; xe=%20.8f nOfBoundStates=%10d\n", toInverseCm(omega0), xe, nOfBoundStates);

        double xmin = 1.0;
        double xmax = 12.;
//        double step = 0.005;
//        int nOfPoints = 1 + (int) ((xmax - xmin) / step);
        System.out.printf("%5s %5s %18s %18s %18s %18s %15s %15s %15s %10s %15s %18s\n", "n", "np", "step", "exact", "calc", "exact-calc", "innerInv", "outerInv",
                "span", "eff.points", "minStep", "step/minStep");
//        System.out.println(toInverseCm(potential.value(xmin)) + " " + toInverseCm(potential.value(xmax)));

        Integrator integrator = IntegratorFactory.getCFMagnus4();

        for (int nOfPoints = 5000; nOfPoints <= 5000; nOfPoints += 1) {
            double step = (xmax - xmin) / (nOfPoints - 1);

            Grid grid = GridFactory.generateUniformGrid(xmin, xmax, nOfPoints);
            SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
//            System.out.printf("minimum step size for all states = %20.8f \n", system.hCriticalAllowed);

            ShootingSolver finder = new ShootingSolver(system, integrator);


//        for (int nOfDesiredNodes = 0; nOfDesiredNodes < nOfBoundStates; nOfDesiredNodes++) {
            for (int nOfDesiredNodes = 100; nOfDesiredNodes < 101; nOfDesiredNodes++) {
                QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes, RefinementStrategy.BISECTION_THEN_BIDIRECTIONAL);
                double exact = omega0 * (nOfDesiredNodes + 0.5) * (1. - xe * (nOfDesiredNodes + 0.5));

                double innerInversionPoint = rmin - Math.log(1. + Math.sqrt(exact / De)) / a;
                double outerInversionPoint = rmin - Math.log(1. - Math.sqrt(exact / De)) / a;
                double span = outerInversionPoint - innerInversionPoint;
                int nEffPoints = (int) (span / step);
                double diff = exact - ek.energy;
                double maxStep = ek.maximumStepSize();
                String msg = (step < maxStep) ? "OK" : "!";
                System.out.printf("%5d %5d %18.8f %18.8f %18.8f %18.10f %15.6f %15.6f %15.6f %10d %18.8f %15.3f %2s\n", nOfDesiredNodes, nOfPoints, step, toInverseCm(exact),
                        toInverseCm(ek.energy), toInverseCm(diff), innerInversionPoint, outerInversionPoint, span, nEffPoints, maxStep, step / maxStep, msg);
            }
            System.out.println("Cache hits   = " + system.getCacheUTilde().nOfCacheHits);
            System.out.println("Cache misses = " + system.getCacheUTilde().nOfCacheMisses);
        }


    }

}
