package lennard_jones;

import org.junit.jupiter.api.Test;
import schrodinger.OutputManager;
import schrodinger.QuantumLevel;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.IntegratorFactory;
import schrodinger.integrator.multi_step.Obrechkoff6;
import schrodinger.integrator.multi_step.Verlet;
import schrodinger.integrator.multi_step.numerovlike.EFN;
import schrodinger.integrator.multi_step.numerovlike.EFNFixedBeta;
import schrodinger.integrator.multi_step.numerovlike.Numerov;
import schrodinger.integrator.multi_step.predcorr.PC6;
import schrodinger.integrator.multi_step.predcorr.PC8i1;
import schrodinger.integrator.multi_step.predcorr.PC8i2;
import schrodinger.integrator.one_step.*;
import schrodinger.integrator.multi_step.Cowell5;
import schrodinger.integrator.multi_step.Cowell6;
import schrodinger.integrator.multi_step.Stormer8;
import schrodinger.integrator.multi_step.Cowell;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialLennardJones;
import schrodinger.SchrodingerSystem;
import schrodinger.solver.RefinementStrategy;
import schrodinger.solver.ShootingSolver;

import java.io.*;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.List;
import java.util.Objects;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static schrodinger.PhysicalConstants.*;

public abstract class LennardJonesAbstractTest {
    private static final Path PROJECT_ROOT = Path.of(System.getProperty("user.dir"));
    private static final Path TEST_SRC_ROOT = PROJECT_ROOT.resolve("test-integration/lennard_jones");

    private static final double wellDepthInverseCm = 4050;
    private static final double wellDepthHartree = wellDepthInverseCm / HARTREE_TO_INVERSE_CM;
    private static final double rMinAng = 1.;
    private static final double rMinBohr = rMinAng / BOHR_TO_ANG;
    private static final PhysicalPotential potential = new PhysicalPotentialLennardJones(rMinBohr, wellDepthHartree);
    private static final double mass = 16.85762920 * UMA_TO_ELECTRON_MASS;

    List<Integrator> integrators;
    private final RefinementStrategy strategy;
    private final boolean isPrintOnlyBad;
    private final Map<String, Double> refEnergiesNew = new HashMap<>();

    void loadReferenceEnergies(String filePath) {

        InputStream is = LennardJonesAbstractTest.class.getResourceAsStream(filePath);
        Objects.requireNonNull(is, "Test resource not found: " + filePath);

        try (is; BufferedReader br = new BufferedReader(new InputStreamReader(is, StandardCharsets.UTF_8))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] parts = line.split("\t");
                String key = parts[0];
                refEnergiesNew.put(key, Double.parseDouble(parts[1]));
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        } catch (NumberFormatException e) {
            System.err.println("Invalid double value in file " + filePath);
        }

    }

    public LennardJonesAbstractTest(RefinementStrategy strategy, boolean isPrintOnlyBad) {
//        this.integrators = IntegratorFactory.getAll();
        this.integrators = List.of(IntegratorFactory.getEFN());
        this.strategy = strategy;
        this.isPrintOnlyBad = isPrintOnlyBad;
    }


    @Test
    void test_1800_points_1_0_to_13_uniform_grid() {
        test_core(1.0, 13, 1800);
    }

    @Test
    void test_1800_points_1_3_to_13_uniform_grid() {
        test_core(1.3, 13, 1800);
    }

    @Test
    void test_1800_points_1_5_to_13_uniform_grid() {
        test_core(1.5, 13, 1800);
    }

    @Test
    void test_1800_points_1_5_to_20_uniform_grid() {
        test_core(1.5, 20, 1800);
    }

    @Test
    void test_1800_points_1_5_to_30_uniform_grid() {
        test_core(1.5, 30, 1800);
    }

    @Test
    void test_1800_points_1_5_to_40_uniform_grid() {
        test_core(1.5, 40, 1800);
    }

    @Test
    void test_1800_points_1_5_to_45_uniform_grid() {
        test_core(1.5, 45, 2000);
    }

    private void test_core(double xmin, double xmax, int nOfPoints) {

        Grid grid = GridFactory.generateUniformGrid(xmin, xmax, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);

        double thresholdAbsInverseCm = 1e-5;
        System.out.printf(String.format("ThresholdAbsInverseCm = %22.3e cm-1\n", thresholdAbsInverseCm));

        int nBad = 0;
        int nGood = 0;
        boolean isFirstRow = true;

        for (Integrator integrator : this.integrators) {
            ShootingSolver finder = new ShootingSolver(system, integrator);

            String className = integrator.getClass().getSimpleName();
            String logFilename = "LennardJones_" + className + "_" + nOfPoints + "_" + this.strategy.toString().toLowerCase() + ".log";
            Path logFile = TEST_SRC_ROOT.resolve(logFilename);
            OutputManager.initCommonOutputFile(logFile.toString());

            System.out.printf("%22s %20s %12s %15s %15s %12s %10s %22s %22s %20s %20s %12s %s\n", "className", "strategy", "nOfPoints", "xmin", "xmax", "nOfDesiredNodes",
                    "goodOrBad", "energy", "energy_ref", "errorAbs", "errorRel", "TotalScans", "<- of which...");

            for (int nOfDesiredNodes = 0; nOfDesiredNodes <= 14; nOfDesiredNodes++) {
                String key = generateKey(className, grid, nOfDesiredNodes);
                Double refEnergy = refEnergiesNew.get(key);
                if (refEnergy == null) {
                    refEnergy = 0.;
                }

                QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes, this.strategy);
                double errorAbs = (refEnergy - toInverseCm(ek.energy));
                double errorRel = errorAbs / refEnergy;

                String goodOrBad;
                if (Math.abs(errorAbs) <= thresholdAbsInverseCm) {
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


                    System.out.printf("%22s %20s %10d %8.2f %8.2f %5d %8s %18.8f %18.8f %18.8f %15.3e %12d %s \n", className, strategy, nOfPoints,
                            xmin, xmax, nOfDesiredNodes,
                            goodOrBad, toInverseCm(ek.energy), refEnergy,
                            errorAbs, errorRel, ek.countTotalScans(), iterInfo);
                }
            }
        }

        double nTotal = (nBad + nGood);
        nTotal = nTotal / 100.;
        System.out.printf("%10s, %10d -- %10.2f%s\n", "nGood: ", nGood, nGood / nTotal, "%");
        System.out.printf("%10s, %10d -- %10.2f%s\n", "nBad", nBad, nBad / nTotal, "%");
        assertEquals(0, nBad);
    }

    private String generateKey(String integratorClassName, Grid grid, int nOfNodes) {
        StringBuilder sb = new StringBuilder();
        sb.append(integratorClassName).append("_");

        String mapping = grid.getMappingStrategy().getClass().getSimpleName();
        if (mapping.isEmpty()) {
            mapping = "uniform";
        }
        sb.append(mapping).append("_");
        sb.append(grid.getNumberOfPoints()).append("_");
        sb.append(grid.getFirstYValue()).append("_");
        sb.append((int) grid.getLastYValue()).append("_");
        sb.append(nOfNodes);
        return sb.toString();
    }

}
