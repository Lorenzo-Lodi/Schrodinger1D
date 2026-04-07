package lennard_jones;

import org.junit.jupiter.api.Test;
import schrodinger.OutputManager;
import schrodinger.QuantumLevel;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.IntegratorFactory;
import schrodinger.integrator.misc.Obrechkoff6;
import schrodinger.integrator.misc.TaylorThreePoints;
import schrodinger.integrator.numerovlike.EFN;
import schrodinger.integrator.numerovlike.EFNFixedBeta;
import schrodinger.integrator.numerovlike.Numerov;
import schrodinger.integrator.predcorr.PredictorCorrector6;
import schrodinger.integrator.predcorr.PredictorCorrector8i1;
import schrodinger.integrator.predcorr.PredictorCorrector8i2;
import schrodinger.integrator.rungekutta.*;
import schrodinger.integrator.stormer.Stormer5;
import schrodinger.integrator.stormer.Stormer6;
import schrodinger.integrator.stormer.Stormer8;
import schrodinger.integrator.stormer.Stormer8i;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialLennardJones;
import schrodinger.potential.SchrodingerSystem;
import schrodinger.solver.RefinementStrategy;
import schrodinger.solver.ShootingSolver;

import java.util.*;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static schrodinger.PhysicalConstants.*;

public class LennardJonesAllTest {

    private static final double wellDepthInverseCm = 4050;
    private static final double wellDepthHartree = wellDepthInverseCm / HARTREE_TO_INVERSE_CM;
    private static final double rMinAng = 1.;
    private static final double rMinBohr = rMinAng / BOHR_TO_ANG;
    private static final PhysicalPotential potential = new PhysicalPotentialLennardJones(rMinBohr, wellDepthHartree, 6);
    private static final double mass = 16.85762920 * UMA_TO_ELECTRON_MASS;
    private static final Map<Integer, Double> refEnergies = new HashMap<>();
    private static final Map<Class<?>, Double> thresholds1800pts = new HashMap<>();

    static {
        // There were obtained with Magnus8, 4000 points and [1.2 - 45.0] uniform grid
        refEnergies.put(0, 371.40006854105536);
        refEnergies.put(1, 1046.35039758025960);
        refEnergies.put(2, 1636.86240202073500);
        refEnergies.put(3, 2147.54058595915200);
        refEnergies.put(4, 2583.19756725932500);
        refEnergies.put(5, 2948.86947887610630);
        refEnergies.put(6, 3249.83161605343780);
        refEnergies.put(7, 3491.61374341683600);
        refEnergies.put(8, 3680.01417388285970);
        refEnergies.put(9, 3821.11132786356030);
        refEnergies.put(10, 3921.27096821614170);
        refEnergies.put(11, 3987.14669051490550);
        refEnergies.put(12, 4025.67057896758570);
        refEnergies.put(13, 4044.03031326710800);
        refEnergies.put(14, 4049.62843497018100);

        thresholds1800pts.put(CFMagnus6e5Opt.class, 3.15E-07);
        thresholds1800pts.put(CFMagnus8.class, 3.15E-07);
        thresholds1800pts.put(PredictorCorrector8i1.class, 3.87E-07);
        thresholds1800pts.put(PredictorCorrector8i2.class, 4.31E-07);
        thresholds1800pts.put(Obrechkoff6.class, 6.70E-07);
        thresholds1800pts.put(CFMagnus6e4.class, 1.36E-06);
        thresholds1800pts.put(Stormer8i.class, 2.62E-06);
        thresholds1800pts.put(PredictorCorrector6.class, 1.31E-05);
        thresholds1800pts.put(CFMagnus4.class, 1.61E-05);
        thresholds1800pts.put(RK45DP.class, 1.93E-05);
        thresholds1800pts.put(Stormer8.class, 4.83E-05);
        thresholds1800pts.put(Stormer5.class, 5.89E-05);
        thresholds1800pts.put(Stormer6.class, 8.00E-05);
        thresholds1800pts.put(EFNFixedBeta.class, 2.33E-04);
        thresholds1800pts.put(EFN.class, 2.33E-04);
        thresholds1800pts.put(Numerov.class, 2.61E-03);
        thresholds1800pts.put(RKN4.class, 3.89E-03);
        thresholds1800pts.put(TaylorThreePoints.class, 1.73E+00);

    }


    @Test
    void testBisectionOnly_1800_points_1_5_to_13_uniform_grid() {
        int nBad = 0;
        int nGood = 0;
        RefinementStrategy strategy = RefinementStrategy.BISECTION_ONLY;
        List<Integrator> integrators = IntegratorFactory.getAll();
        for (Integrator integrator : integrators) {
            int nOfPoints = 1800;
            for (int nOfDesiredNodes = 0; nOfDesiredNodes <= 14; nOfDesiredNodes++) {
                String className = integrator.getClass().getSimpleName().replaceAll("Test", "");
                OutputManager.initCommonOutputFile("LennardJones_" + className + "_" +
                        nOfPoints + "_" + strategy.toString().toLowerCase() + ".log");
                Grid grid = GridFactory.generateUniformGrid(1.5d, 13., nOfPoints);
                OutputManager.write("Grid info");
                OutputManager.write(String.format("%10s %22s %22s %22s", "i", "r", "y", "V"));
//                for (int i = 0; i < grid.getNumberOfPoints(); i++) {
//                    double v = potential.value(grid.rAtGridPoint(i));
//                    OutputManager.write(String.format("%10d %22.8f %22.8f %25.6f",
//                            i, grid.rAtGridPoint(i), grid.yAtGridPoint(i), toInverseCm(v)));
//                }
                SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
                ShootingSolver finder = new ShootingSolver(system, integrator);
                QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes, strategy);
                double error = (refEnergies.get(nOfDesiredNodes) - toInverseCm(ek.energy)) / refEnergies.get(nOfDesiredNodes);
                double threshold = thresholds1800pts.get(integrator.getClass());
                if(nOfDesiredNodes == 14) {
                    threshold = 0.07 ; // Special rule for the highest level, as it is not completely converged because of grid
                }
                String goodOrBad;
                if (Math.abs(error) <= threshold) {
                    goodOrBad = "Good";
                    nGood++;
                } else {
                    goodOrBad = "Bad";
                    nBad++;
                }
                System.out.printf("%22s %15s %8d %8d %10s %20.4f %20.4f %20.8e \n", className, strategy, nOfPoints, nOfDesiredNodes,
                        goodOrBad, toInverseCm(ek.energy), refEnergies.get(nOfDesiredNodes), error);
            }
        }

        double nTotal = (nBad + nGood);
        nTotal = nTotal / 1000.;
        System.out.printf("%10s, %10d -- %10.2f%s\n", "nGood: ", nGood, nGood / nTotal, "%");
        System.out.printf("%10s, %10d -- %10.2f%s\n", "nBad", nBad, nBad / nTotal, "%");
        assertEquals(0, nBad);
    }


    @Test
    void testBisectionOnly_1800_points_1_5_to_20_uniform_grid() {
        int nBad = 0;
        int nGood = 0;
        RefinementStrategy strategy = RefinementStrategy.BISECTION_ONLY;
        List<Integrator> integrators = IntegratorFactory.getAll();
//        List<Integrator> integrators = List.of(IntegratorFactory.getNumerov());
        for (Integrator integrator : integrators) {
            int nOfPoints = 1800;
            for (int nOfDesiredNodes = 0; nOfDesiredNodes <= 14; nOfDesiredNodes++) {
                String className = integrator.getClass().getSimpleName().replaceAll("Test", "");
                OutputManager.initCommonOutputFile("LennardJones_" + className + "_" +
                        nOfPoints + "_" + strategy.toString().toLowerCase() + ".log");
                Grid grid = GridFactory.generateUniformGrid(1.5d, 20., nOfPoints);
                OutputManager.write("Grid info");
                OutputManager.write(String.format("%10s %22s %22s %22s", "i", "r", "y", "V"));
//                for (int i = 0; i < grid.getNumberOfPoints(); i++) {
//                    double v = potential.value(grid.rAtGridPoint(i));
//                    OutputManager.write(String.format("%10d %22.8f %22.8f %25.6f",
//                            i, grid.rAtGridPoint(i), grid.yAtGridPoint(i), toInverseCm(v)));
//                }
                SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
                ShootingSolver finder = new ShootingSolver(system, integrator);
                QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes, strategy);
                double error = (refEnergies.get(nOfDesiredNodes) - toInverseCm(ek.energy)) / refEnergies.get(nOfDesiredNodes);
                double threshold = thresholds1800pts.get(integrator.getClass());
                if(nOfDesiredNodes == 14) {
                    threshold = 0.001 ; // Special rule for the highest level, as it is not completely converged because of grid
                }
                String goodOrBad;
                if (Math.abs(error) <= threshold) {
                    goodOrBad = "Good";
                    nGood++;
                } else {
                    goodOrBad = "Bad";
                    nBad++;
                }
                System.out.printf("%22s %15s %8d %8d %10s %20.4f %20.4f %20.8e \n", className, strategy, nOfPoints, nOfDesiredNodes,
                        goodOrBad, toInverseCm(ek.energy), refEnergies.get(nOfDesiredNodes), error);
            }
        }

        double nTotal = (nBad + nGood);
        nTotal = nTotal / 1000.;
        System.out.printf("%10s, %10d -- %10.2f%s\n", "nGood: ", nGood, nGood / nTotal, "%");
        System.out.printf("%10s, %10d -- %10.2f%s\n", "nBad", nBad, nBad / nTotal, "%");
        assertEquals(0, nBad);
    }


}
