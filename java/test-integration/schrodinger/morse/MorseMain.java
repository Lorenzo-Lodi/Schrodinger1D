package schrodinger.morse;

import schrodinger.OutputManager;
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

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static schrodinger.PhysicalConstants.*;
import static schrodinger.Utils.stepMsg;

public class MorseMain {
    private static final Path PROJECT_ROOT = Path.of(System.getProperty("user.dir"));
    private static final Path TEST_SRC_ROOT = PROJECT_ROOT.resolve("test-integration/schrodinger/morse");
    private final static Map<String, Double> refEnergies = new HashMap<>();

    public static void main(String[] args) {
        morse3();
    }


    private static void morse_core() {

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
        double step = 0.005;
        int nOfPoints = 1 + (int) ((xmax - xmin) / step);
        System.out.printf("%15s %5s %5s %18s %18s %18s %18s %15s %15s %15s %10s %15s %18s\n", "ClassName", "n", "np", "step", "exact", "calc", "exact-calc", "innerInv", "outerInv",
                "span", "eff.points", "minStep", "step/minStep");
        System.out.println(toInverseCm(potential.value(xmin)) + " " + toInverseCm(potential.value(xmax)));

        RefinementStrategy strategy = RefinementStrategy.BISECTION_THEN_BIDIRECTIONAL;

        for (Integrator integrator : List.of(IntegratorFactory.getCFMagnus4())) {
            int totalScans = 0;
            String className = integrator.getClass().getSimpleName();
            long t0 = System.nanoTime();

            String logFilename = "Morse_" + className + "_" + xmin + "_" + xmax + "_" + nOfPoints + "_"
                    + strategy.toString().toLowerCase() + ".log";
            Path logFile = TEST_SRC_ROOT.resolve("outputs/" + logFilename);
            OutputManager.initCommonOutputFile(logFile.toString());

            Grid grid = GridFactory.generateUniformGrid(xmin, xmax, nOfPoints);
            SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
            system.setCachingUTilde(true);
            int nOffsets = integrator.getFractionalOffsets().length;
            System.out.printf("Maximum step size for ALL states to dissociation = %20.8f \n", system.maxStepSizeAllowedRegion / nOffsets);
            ShootingSolver finder = new ShootingSolver(system, integrator);

            for (int nOfDesiredNodes = 0; nOfDesiredNodes < 101; nOfDesiredNodes++) {
                QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes, strategy);
                double exact = omega0 * (nOfDesiredNodes + 0.5) * (1. - xe * (nOfDesiredNodes + 0.5));

                double innerInversionPoint = rmin - Math.log(1. + Math.sqrt(exact / De)) / a;
                double outerInversionPoint = rmin - Math.log(1. - Math.sqrt(exact / De)) / a;
                double span = outerInversionPoint - innerInversionPoint;
                int nEffPoints = (int) (span / step);
                double diff = exact - ek.energy;
                double maxStep = ek.maximumStepSize() * nOffsets;
                totalScans += ek.countTotalScans();
                String msg = stepMsg(step / maxStep);
                System.out.printf("%15s %5d %5d %18.8f %18.8f %18.8f %18.10f %15.6f %15.6f %15.6f %10d %18.8f %15.3f %2s %6d\n", className, nOfDesiredNodes, nOfPoints, step, toInverseCm(exact),
                        toInverseCm(ek.energy), toInverseCm(diff), innerInversionPoint, outerInversionPoint, span, nEffPoints, maxStep, step / maxStep, msg, ek.countTotalScans());
            }
            int hits = system.getCacheUTilde().nOfCacheHits;
            int misses = system.getCacheUTilde().nOfCacheMisses;
            double totCache = hits + misses;
            double hitRate = (totCache == 0) ? 0.0 : 100. * hits / totCache;
            double missRate = (totCache == 0) ? 0.0 : 100. * misses / totCache;
//            System.out.printf(
//                    "Cache hits = %12d (%7.3f %%), misses = %12d (%7.3f %%)%n",
//                    hits, hitRate, misses, missRate
//            );
            long elapsedNanos = System.nanoTime() - t0;
            double elapsedSeconds = elapsedNanos * 1e-9;
            System.out.printf("%15s wall-time seconds = %25.6f ; Cache hit rate = %7.3f %%; Total scans = %12d \n", className, elapsedSeconds,
                    hitRate, totalScans);
        }
    }

    private static void morse2() {

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
        double step = 0.005;
        int nOfPoints = 1 + (int) ((xmax - xmin) / step);
        System.out.printf("%15s %5s %5s %18s %18s %18s %18s %15s %15s %15s %10s %15s %18s\n", "ClassName", "n", "np", "step", "LEVEL", "calc", "LEVEL-calc", "innerInv", "outerInv",
                "span", "eff.points", "minStep", "step/minStep");
        System.out.println(toInverseCm(potential.value(xmin)) + " " + toInverseCm(potential.value(xmax)));

        RefinementStrategy strategy = RefinementStrategy.BISECTION_THEN_BIDIRECTIONAL;

        for (Integrator integrator : List.of(IntegratorFactory.getCFMagnus4())) {
            int totalScans = 0;
            String className = integrator.getClass().getSimpleName();
            long t0 = System.nanoTime();

            String logFilename = "Morse_" + className + "_" + xmin + "_" + xmax + "_" + nOfPoints + "_"
                    + strategy.toString().toLowerCase() + ".log";
            Path logFile = TEST_SRC_ROOT.resolve("outputs/" + logFilename);
            OutputManager.initCommonOutputFile(logFile.toString());

            Grid grid = GridFactory.generateUniformGrid(xmin, xmax, nOfPoints);
            SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
            system.setCachingUTilde(true);
            int nOffsets = integrator.getFractionalOffsets().length;
            System.out.printf("Maximum step size for ALL states to dissociation = %20.8f \n", system.maxStepSizeAllowedRegion / nOffsets);
            ShootingSolver finder = new ShootingSolver(system, integrator);

            List<QuantumLevel> levels = finder.findEigenvaluesUpTo(100, strategy);

            for (int nOfDesiredNodes = 0; nOfDesiredNodes < 101; nOfDesiredNodes++) {
                QuantumLevel ek = levels.get(nOfDesiredNodes);
                double exact = omega0 * (nOfDesiredNodes + 0.5) * (1. - xe * (nOfDesiredNodes + 0.5));

                double innerInversionPoint = rmin - Math.log(1. + Math.sqrt(exact / De)) / a;
                double outerInversionPoint = rmin - Math.log(1. - Math.sqrt(exact / De)) / a;
                double span = outerInversionPoint - innerInversionPoint;
                int nEffPoints = (int) (span / step);
                double diff = exact - ek.energy;
                double maxStep = ek.maximumStepSize() * nOffsets;
                totalScans += ek.countTotalScans();
                String msg = stepMsg(step / maxStep);
                System.out.printf("%15s %5d %5d %18.8f %18.8f %18.8f %18.10f %15.6f %15.6f %15.6f %10d %18.8f %15.3f %2s %6d\n", className, nOfDesiredNodes, nOfPoints, step, toInverseCm(exact),
                        toInverseCm(ek.energy), toInverseCm(diff), innerInversionPoint, outerInversionPoint, span, nEffPoints, maxStep, step / maxStep, msg, ek.countTotalScans());
            }


            int hits = system.getCacheUTilde().nOfCacheHits;
            int misses = system.getCacheUTilde().nOfCacheMisses;
            double totCache = hits + misses;
            double hitRate = (totCache == 0) ? 0.0 : 100. * hits / totCache;
            double missRate = (totCache == 0) ? 0.0 : 100. * misses / totCache;
//            System.out.printf(
//                    "Cache hits = %12d (%7.3f %%), misses = %12d (%7.3f %%)%n",
//                    hits, hitRate, misses, missRate
//            );
            long elapsedNanos = System.nanoTime() - t0;
            double elapsedSeconds = elapsedNanos * 1e-9;
            System.out.printf("%15s wall-time seconds = %25.6f ; Cache hit rate = %7.3f %%; Total scans = %12d \n", className, elapsedSeconds,
                    hitRate, totalScans);
        }
    }

    private static void morse3() {

        String filePath = "./resources/level16_reference_values.txt";
        Path inputFile = TEST_SRC_ROOT.resolve(filePath);

        try (InputStream is = Files.newInputStream(inputFile);
             BufferedReader br = new BufferedReader(new InputStreamReader(is, StandardCharsets.UTF_8))) {
            String line;
            br.readLine(); // Skip line
            br.readLine(); // Skip line
            while ((line = br.readLine()) != null) {
                line = line.replace("  ", " ");
                line = line.replace("  ", " ");
                line = line.replace("  ", " ");
                line = line.replace("  ", " ");
                String[] parts = line.split(" ");
                String key = parts[1] + "_" + parts[3];
                refEnergies.put(key, Double.parseDouble(parts[1]));
            }
        } catch (IOException e) {
            throw new UncheckedIOException("Test input file not found: " + inputFile, e);
        } catch (NumberFormatException e) {
            System.err.println("Invalid double value in file " + filePath);
        }
        System.out.println("Loaded " + refEnergies.size() + " reference energies from " + inputFile);



        double a = 1.5;
        double rmin = 2.;
        double De = toHartree(50000.);

        PhysicalPotential potential = new PhysicalPotentialMorse(rmin, a, De);
        double mass = 100. * UMA_TO_ELECTRON_MASS;

        double omega0 = a * Math.sqrt(2.0 * De / mass);
        double xe = omega0 / (4. * De);
        double A = 1. / xe;
        int nOfBoundStates = (int) ((A + 1.) * 0.5);

        double xmin = 1.0;
        double xmax = 12.;
        double step = 0.001;
        int nOfPoints = 1 + (int) ((xmax - xmin) / step);
        System.out.printf("%15s %5s %5s %5s %18s %18s %18s %18s %15s %15s %15s %10s %15s %18s\n", "ClassName", "jrot", "v", "np", "step", "exact", "calc", "exact-calc", "innerInv", "outerInv",
                "span", "eff.points", "minStep", "step/minStep");

        RefinementStrategy strategy = RefinementStrategy.BISECTION_THEN_BIDIRECTIONAL;

//        level16_reference_values.txt

        Integrator integrator = IntegratorFactory.getNumerov();
        int totalScans = 0;
        String className = integrator.getClass().getSimpleName();
        long t0 = System.nanoTime();

        String logFilename = "Morse_" + className + "_" + xmin + "_" + xmax + "_" + nOfPoints + "_"
                + strategy.toString().toLowerCase() + ".log";
        Path logFile = TEST_SRC_ROOT.resolve("outputs/" + logFilename);
        OutputManager.initCommonOutputFile(logFile.toString());

        Grid grid = GridFactory.generateUniformGrid(xmin, xmax, nOfPoints);
        SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
        system.setCachingUTilde(true);
        int nOffsets = integrator.getFractionalOffsets().length;
//            System.out.printf("Maximum step size for ALL states to dissociation = %20.8f \n", system.maxStepSizeAllowedRegion / nOffsets);
        ShootingSolver finder = new ShootingSolver(system, integrator);

        int vMax = 0;
        int Jmax = 10;
        List<List<QuantumLevel>> levels = finder.findEigenvaluesForJ(vMax, Jmax);
        for (int jrot = 0; jrot < levels.size(); jrot++) {
            List<QuantumLevel> levelsPerJ = levels.get(jrot);
            for (int v = 0; v < levelsPerJ.size(); v++) {
                QuantumLevel ek = levelsPerJ.get(v);
                double exact = 0.;
                double innerInversionPoint = rmin - Math.log(1. + Math.sqrt(exact / De)) / a;
                double outerInversionPoint = rmin - Math.log(1. - Math.sqrt(exact / De)) / a;
                double span = outerInversionPoint - innerInversionPoint;
                int nEffPoints = (int) (span / step);
                double diff = exact - ek.energy;
                double maxStep = ek.maximumStepSize() * nOffsets;
                totalScans += ek.countTotalScans();
                String msg = stepMsg(step / maxStep);
                System.out.printf("%15s %5d %5d %5d %18.8f %18.8f %18.8f %18.10f %15.6f %15.6f %15.6f %10d %18.8f %15.3f %2s %6d\n", className, v, jrot, nOfPoints, step, toInverseCm(exact),
                        toInverseCm(ek.energy), toInverseCm(diff), innerInversionPoint, outerInversionPoint, span, nEffPoints, maxStep, step / maxStep, msg, ek.countTotalScans());
            }

            int hits = system.getCacheUTilde().nOfCacheHits;
            int misses = system.getCacheUTilde().nOfCacheMisses;
            double totCache = hits + misses;
            double hitRate = (totCache == 0) ? 0.0 : 100. * hits / totCache;
            double missRate = (totCache == 0) ? 0.0 : 100. * misses / totCache;
//            System.out.printf(
//                    "Cache hits = %12d (%7.3f %%), misses = %12d (%7.3f %%)%n",
//                    hits, hitRate, misses, missRate
//            );
            long elapsedNanos = System.nanoTime() - t0;
            double elapsedSeconds = elapsedNanos * 1e-9;
            System.out.printf("%15s wall-time seconds = %25.6f ; Cache hit rate = %7.3f %%; Total scans = %12d \n", className, elapsedSeconds,
                    hitRate, totalScans);
        }
    }

}
