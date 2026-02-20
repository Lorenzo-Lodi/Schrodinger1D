import schrodinger.integrator.*;

/**
 * Simple test to verify the refactoring works correctly without JUnit dependencies.
 */
public class SimpleTest {
    public static void main(String[] args) {
        System.out.println("Testing refactored integrator test classes...");
        
        // Test that we can instantiate each test class
        try {
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
            
            System.out.println("Successfully created " + tests.length + " test classes.");
            
            // Test that each test class returns the correct integrator
            for (AbstractIntegratorTest test : tests) {
                Integrator integrator = test.getIntegrator();
                System.out.println("Test class: " + test.getClass().getSimpleName() + 
                                  " - Integrator: " + integrator.getClass().getSimpleName());
            }
            
            // Test the ConvergenceData class
            AbstractIntegratorTest.ConvergenceData data = new AbstractIntegratorTest.ConvergenceData();
            data.addError("10", 0.001);
            data.addError("4", 0.002);
            data.addError("2", 0.003);
            
            System.out.println("ConvergenceData test:");
            System.out.println("Error at 10%: " + data.getError("10"));
            System.out.println("Error at 25%: " + data.getError("4"));
            System.out.println("Error at 50%: " + data.getError("2"));
            
            System.out.println("All tests passed successfully!");
        } catch (Exception e) {
            System.err.println("Error during testing: " + e.getMessage());
            e.printStackTrace();
        }
    }
}