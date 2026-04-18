package lennard_jones;

import schrodinger.OutputManager;
import schrodinger.QuantumLevel;
import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.integrator.Integrator;
import schrodinger.integrator.IntegratorFactory;
import schrodinger.potential.PhysicalPotential;
import schrodinger.potential.PhysicalPotentialLennardJones;
import schrodinger.potential.SchrodingerSystem;
import schrodinger.solver.RefinementStrategy;
import schrodinger.solver.ShootingSolver;

import java.util.HashMap;
import java.util.Map;

import static schrodinger.PhysicalConstants.*;

public class LennardJonesMain {
    private static final double wellDepthInverseCm = 4050;
    private static final double wellDepthHartree = wellDepthInverseCm / HARTREE_TO_INVERSE_CM;
    private static final double rMinAng = 1.;
    private static final double rMinBohr = rMinAng / BOHR_TO_ANG;
    private static final PhysicalPotential potential = new PhysicalPotentialLennardJones(rMinBohr, wellDepthHartree, 6);
    private static final double mass = 16.85762920 * UMA_TO_ELECTRON_MASS;
    private static final Map<Integer, Double> refEnergies = new HashMap<>();

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
    }

    public static void main(String[] args) {
        double xmin = 1.4;
        int nOfPoints = 1800;
        Integrator integrator = IntegratorFactory.getStormer8();

        int nOfDesiredNodes = 9;
        RefinementStrategy strategy = RefinementStrategy.BISECTION_ONLY;
        System.out.printf("%20s %20s %25s\n", "xmax", "ptDensity", "errorAbs");
        for (int i = 20; i <= 60; i++) {
            OutputManager.initCommonOutputFile("LennardJones_" + "Stormer8" + "_" +
                    nOfPoints + "_" + strategy.toString().toLowerCase() + ".log");
            double xmax = i;
            double stepSize = (xmax - xmin) / (nOfPoints - 1);
            double ptDensity = 1. / stepSize;
            Grid grid = GridFactory.generateUniformGrid(xmin, xmax, nOfPoints);
            SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
            ShootingSolver finder = new ShootingSolver(system, integrator);
            QuantumLevel ek = finder.findEigenvalue(nOfDesiredNodes, strategy);
            double errorAbs = (refEnergies.get(nOfDesiredNodes) - toInverseCm(ek.energy));
            System.out.printf("%20.4f %20.4f %25.8f\n", xmax, ptDensity, errorAbs);
        }


    }

}
