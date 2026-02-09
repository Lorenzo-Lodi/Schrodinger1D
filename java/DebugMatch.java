import schrodinger.grid.Grid;
import schrodinger.grid.GridFactory;
import schrodinger.potential.PhysicalPotentialHarmonic;
import schrodinger.potential.SchrodingerSystem;

public class DebugMatch {
    public static void main(String[] args) {
        PhysicalPotentialHarmonic potential = new PhysicalPotentialHarmonic(0, 1);
        double mass = 2.0;
        double energy = 0.676; // approximate converged energy
        
        System.out.println("nPoints   h         matchIdx   matchY");
        for (int nPoints = 20; nPoints <= 550; nPoints += 50) {
            Grid grid = GridFactory.generateUniformGrid(-5.0, 5.0, nPoints);
            SchrodingerSystem system = new SchrodingerSystem(potential, mass, grid);
            system.setEnergy(energy);
            
            int matchIdx = findMatchingIndex(system, energy);
            double matchY = grid.getYValue(matchIdx);
            System.out.println(nPoints + "       " + String.format("%.4f", grid.getStepSizeYCoordinate()) + 
                             "     " + matchIdx + "       " + String.format("%.4f", matchY));
        }
    }
    
    private static int findMatchingIndex(SchrodingerSystem system, double energy) {
        var grid = system.getGrid();
        for (int i = grid.getNumberOfPoints() - 3; i >= 2; i--) {
            if (system.UTilde(grid.getYValue(i)) <= energy) {
                return i;
            }
        }
        return grid.getNumberOfPoints() / 2;
    }
}
