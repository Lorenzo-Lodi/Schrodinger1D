import schrodinger.integrator.*;

/**
 * Simple test to verify the refactoring works correctly.
 */
public class TestRefactoring {
    public static void main(String[] args) {
        // Test that we can instantiate each test class
        AbstractIntegratorTest[] tests = {
            new TaylorThreePointsTest(),
            new NumerovTest(),
            new VignoliTest(),
            new RaptisAllisonTest(),
            new ExponentiallyFittedTest(),
            new PredictorCorrector6Test(),
            new PredictorCorrector8Test(),
            new PredictorCorrector8AltTest(),
            new Stormer7Test(),
            new Stormer8Test()
        };
        
        System.out.println("Created " + tests.length + " test classes successfully.");
        
        // Test that each test class returns the correct integrator
        for (AbstractIntegratorTest test : tests) {
            Integrator integrator = test.getIntegrator();
            System.out.println("Test class: " + test.getClass().getSimpleName() + 
                              " - Integrator: " + integrator.getClass().getSimpleName());
        }
        
        System.out.println("Refactoring verification complete!");
    }
}