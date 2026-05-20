package lennard_jones;

import org.junit.jupiter.api.Test;
import schrodinger.OutputManager;
import schrodinger.QuantumLevel;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.IntegratorFactory;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialLennardJones;
import schrodinger.SchrodingerSystem;
import schrodinger.solver.RefinementStrategy;
import schrodinger.solver.ShootingSolver;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import java.util.List;

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
    private final Map<String, Double> refEnergies = new HashMap<>();

    void loadReferenceEnergies(String filePath) {
        Path inputFile = TEST_SRC_ROOT.resolve(filePath);

        try (InputStream is = Files.newInputStream(inputFile);
             BufferedReader br = new BufferedReader(new InputStreamReader(is, StandardCharsets.UTF_8))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] parts = line.split("\t");
                String key = parts[0];
                refEnergies.put(key, Double.parseDouble(parts[1]));
            }
        } catch (IOException e) {
            throw new UncheckedIOException("Test input file not found: " + inputFile, e);
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

        String iterDescription = "";

        for (Integrator integrator : this.integrators) {
            ShootingSolver finder = new ShootingSolver(system, integrator);

            String className = integrator.getClass().getSimpleName();
            String logFilename = "LennardJones_" + className + "_" + nOfPoints + "_" + this.strategy.toString().toLowerCase() + ".log";
            Path logFile = TEST_SRC_ROOT.resolve("outputs/" + logFilename);
            OutputManager.initCommonOutputFile(logFile.toString());

            System.out.printf("%22s %30s %10s %8s %8s %8s %5s %18s %18s %18s %15s %12s %17s | %8s %8s %8s\n",
                    "className", "strategy", "nOfPoints", "xmin", "xmax", "step",
                    "nodes", "energy", "energy_ref", "errorAbs", "errorRel", "check", " SCANS...",
                    "LftInv", "RghtInv", "minStep");


            for (int nOfDesiredNodes = 0; nOfDesiredNodes <= 14; nOfDesiredNodes++) {
                String key = generateKey(className, grid, nOfDesiredNodes);
                Double refEnergy = refEnergies.get(key);
                if (refEnergy == null) {
                    refEnergy = 0.;
                }

                QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes, this.strategy);
                double errorAbs = (refEnergy - toInverseCm(ek.energy));
                double errorRel = errorAbs / refEnergy;

                double leftInversionPoint = system.findLeftmostInversionGridpointY(toHartree(refEnergy));
                double rightInversionPoint = system.findRightmostInversionGridpointY(toHartree(refEnergy));
                double allowedRange = rightInversionPoint - leftInversionPoint;

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
                    iterDescription = "";
                    for (int j = 0; j < ek.convergenceInfo.size(); j++) {
                        QuantumLevel.ConvergenceInfo info = ek.convergenceInfo.get(j);
                        iterInfo += String.format("%4d", info.iterations);
                        iterDescription += String.format("%3d %s\n", j, info.convergengeStage);
                    }

                    System.out.printf("%22s %30s %10d %8.2f %8.2f %8.4f %5d %18.8f %18.8f %18.8f %15.3e %12s %4d %s | %8.2f %8.2f %8.4f\n",
                            className, strategy, nOfPoints,
                            xmin, xmax, grid.getStepSizeYCoordinate(),
                            nOfDesiredNodes,
                            toInverseCm(ek.energy), refEnergy,
                            errorAbs, errorRel, goodOrBad, ek.countTotalScans(), iterInfo,
                            leftInversionPoint, rightInversionPoint,
                            ek.minimumStepSize());
                }
            }
        }
        System.out.println("SCAN iterations description:");
        System.out.printf(iterDescription + "\n");

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
