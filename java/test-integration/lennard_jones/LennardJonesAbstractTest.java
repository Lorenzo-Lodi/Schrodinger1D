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

import java.util.HashMap;
import java.util.Map;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static schrodinger.PhysicalConstants.*;

public abstract class LennardJonesAbstractTest {

    private static final double wellDepthInverseCm = 4050;
    private static final double wellDepthHartree = wellDepthInverseCm / HARTREE_TO_INVERSE_CM;
    private static final double rMinAng = 1.;
    private static final double rMinBohr = rMinAng / BOHR_TO_ANG;
    private static final PhysicalPotential potential = new PhysicalPotentialLennardJones(rMinBohr, wellDepthHartree, 6);
    private static final double mass = 16.85762920 * UMA_TO_ELECTRON_MASS;
    private static final Map<Integer, Double> refEnergies = new HashMap<>();
    private static final Map<Class<?>, Double> thresholds1800pts = new HashMap<>();

    List<Integrator> integrators;
    private final RefinementStrategy strategy;
    private final boolean isPrintOnlyBad;

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

        // Thresholds for 1800 points, grid 1.5-13; thresholds for all 14 states (worse-case), ABSOLUTE ERRORS
        thresholds1800pts.put(CFMagnus6e5Opt.class, 3.16E-07);
        thresholds1800pts.put(CFMagnus8.class, 3.16E-07);
        thresholds1800pts.put(PredictorCorrector8i1.class, 4.6E-07);
        thresholds1800pts.put(PredictorCorrector8i2.class, 6.5E-07);
        thresholds1800pts.put(Obrechkoff6.class, 1.10E-06);
        thresholds1800pts.put(CFMagnus6e4.class, 3.6E-06);
        thresholds1800pts.put(Stormer8i.class, 3.2E-06);
        thresholds1800pts.put(PredictorCorrector6.class, 1.32E-05);
        thresholds1800pts.put(CFMagnus4.class, 1.61E-05);
        thresholds1800pts.put(RK45DP.class, 6.0E-05);
        thresholds1800pts.put(Stormer8.class, 6.63E-05);
        thresholds1800pts.put(Stormer5.class, 1.7E-04);
        thresholds1800pts.put(Stormer6.class, 8.01E-05);
        thresholds1800pts.put(EFNFixedBeta.class, 2.33E-04);
        thresholds1800pts.put(EFN.class, 2.33E-04);
        thresholds1800pts.put(Numerov.class, 2.61E-03);
        thresholds1800pts.put(RKN4.class, 3.89E-03);
        thresholds1800pts.put(TaylorThreePoints.class, 1.74E+00);

    }

    public LennardJonesAbstractTest(RefinementStrategy strategy, boolean isPrintOnlyBad) {
        this.integrators = IntegratorFactory.getAll();
        this.strategy = strategy;
        this.isPrintOnlyBad = isPrintOnlyBad;
    }


    @Test
    void test_1800_points_1_0_to_13_uniform_grid() {
        test_core(1.0, 13, 1800, 0.07);
    }

    @Test
    void test_1800_points_1_3_to_13_uniform_grid() {
        test_core(1.3, 13, 1800, 0.07);
    }

    @Test
    void test_1800_points_1_5_to_13_uniform_grid() {
        test_core(1.5, 13, 1800, 0.07);
    }

    @Test
    void test_1800_points_1_5_to_20_uniform_grid() {
        test_core(1.5, 20, 1800, 0.06);
    }

    @Test
    void test_1800_points_1_5_to_30_uniform_grid() {
        test_core(1.5, 30, 1800, 0.14);
    }

    @Test
    void test_1800_points_1_5_to_40_uniform_grid() {
        test_core(1.5, 40, 1800, 0.28);
    }

    @Test
    void test_1800_points_1_5_to_45_uniform_grid() {
        test_core(1.5, 45, 1800, 0.40);
    }


    private void test_core(double xmin, double xmax, int nOfPoints, double threshold14thState) {
        int nBad = 0;
        int nGood = 0;

        boolean isFirstRow = true;

        for (Integrator integrator : this.integrators) {

            if(xmax >= 45. && integrator instanceof Stormer8) {
                continue; //For now skip Stormer 8 for very long bond lengths as it goes crazy
            }

            String className = integrator.getClass().getSimpleName().replaceAll("Test", "");
            OutputManager.initCommonOutputFile("LennardJones_" + className + "_" +
                    nOfPoints + "_" + this.strategy.toString().toLowerCase() + ".log");
            Grid grid = GridFactory.generateUniformGrid(xmin, xmax, nOfPoints);
//                OutputManager.write("Grid info");
//                OutputManager.write(String.format("%10s %22s %22s %22s", "i", "r", "y", "V"));
//                for (int i = 0; i < grid.getNumberOfPoints(); i++) {
//                    double v = potential.value(grid.rAtGridPoint(i));
//                    OutputManager.write(String.format("%10d %22.8f %22.8f %25.6f",
//                            i, grid.rAtGridPoint(i), grid.yAtGridPoint(i), toInverseCm(v)));
//                }
            SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
            ShootingSolver finder = new ShootingSolver(system, integrator);

            double thresholdAbs = thresholds1800pts.get(integrator.getClass());
            thresholdAbs = 1.25 * thresholdAbs * Math.pow((xmax - xmin) / (13. - 1.5), integrator.globalConvergenceOrder());
            System.out.printf(String.format("For class %22s thresholdAbs (max error for all states apart nNodes=14) is currently set to: %22.3e cm-1\n", className, thresholdAbs));
            System.out.printf("%22s %20s %12s %15s %15s %20s %10s %22s %22s %20s %20s %12s %s\n", "className", "strategy", "nOfPoints", "xmin", "xmax", "nOfDesiredNodes",
                    "goodOrBad", "energy", "energy_ref", "errorAbs", "errorRel", "TotalScans", "<- of which...");

            for (int nOfDesiredNodes = 0; nOfDesiredNodes <= 14; nOfDesiredNodes++) {
                QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes, this.strategy);
                double errorAbs = (refEnergies.get(nOfDesiredNodes) - toInverseCm(ek.energy));
                double errorRel = errorAbs / refEnergies.get(nOfDesiredNodes);
                // Special rule for the highest level, as it is not completely converged because of grid
                if (nOfDesiredNodes == 14) {
                    thresholdAbs = threshold14thState;
                }
                String goodOrBad;
                if (Math.abs(errorAbs) <= thresholdAbs) {
                    goodOrBad = "Good";
                    nGood++;
                } else {
                    goodOrBad = "Bad";
                    nBad++;
                }

                if (!(isPrintOnlyBad && goodOrBad.equals("Good"))) {
                    String iterInfo = "";
                    String iterDescription = "";
                    for (int j = 0; j < ek.convergenceInfo.size(); j++) {
                        QuantumLevel.ConvergenceInfo info = ek.convergenceInfo.get(j);
                        iterInfo += String.format("%4d", info.iterations);
                        iterDescription += String.format("%3d %s\n", j, info.convergengeStage);
                    }
                    if (isFirstRow) {
                        System.out.printf(iterDescription);
                        isFirstRow = false;
                    }

                    System.out.printf("%22s %20s %12d %15.8e %15.8e %20d %10s %22.8f %22.8f %20.8f %20.8e %12d %s\n", className, strategy, nOfPoints,
                            xmin, xmax, nOfDesiredNodes,
                            goodOrBad, toInverseCm(ek.energy), refEnergies.get(nOfDesiredNodes), errorAbs, errorRel, ek.countTotalScans(), iterInfo);
                }
            }
        }

        double nTotal = (nBad + nGood);
        nTotal = nTotal / 100.;
        System.out.printf("%10s, %10d -- %10.2f%s\n", "nGood: ", nGood, nGood / nTotal, "%");
        System.out.printf("%10s, %10d -- %10.2f%s\n", "nBad", nBad, nBad / nTotal, "%");
        assertEquals(0, nBad);
    }


}
