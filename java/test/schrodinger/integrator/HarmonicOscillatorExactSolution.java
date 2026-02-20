package schrodinger.integrator;

/**
 * Provides exact analytical solutions for the quantum harmonic oscillator.
 * Supports ground state and any excited state (n >= 0).
 * 
 * The potential is V(r) = alpha * (r - r0)^2
 * The exact wavefunctions are:
 * ψ_n(x) = N_n * H_n(x) * exp(-x^2/2)
 * 
 * where:
 * - N_n = (1 / sqrt(2^n * n! * sqrt(pi))) is the normalization constant
 * - H_n(x) are Hermite polynomials
 * - x = sqrt(alpha) * (r - r0) is the reduced coordinate
 * 
 * Energy levels: E_n = (n + 1/2) * omega (in appropriate units)
 */
public class HarmonicOscillatorExactSolution {
    
    private final double r0;
    private final double alpha;
    private final int quantumNumber;
    private final boolean normalized;
    
    /**
     * Creates an exact solution for the n-th quantum state of a harmonic oscillator.
     * By default, returns unnormalized wavefunction for backward compatibility.
     * 
     * @param r0 The center of the harmonic potential
     * @param alpha The spring constant parameter (V = alpha * (r - r0)^2)
     * @param quantumNumber The quantum number n (0 = ground state, 1 = first excited, etc.)
     */
    public HarmonicOscillatorExactSolution(double r0, double alpha, int quantumNumber) {
        this(r0, alpha, quantumNumber, false); // Default to unnormalized for backward compatibility
    }
    
    /**
     * Creates an exact solution for the n-th quantum state of a harmonic oscillator.
     * 
     * @param r0 The center of the harmonic potential
     * @param alpha The spring constant parameter (V = alpha * (r - r0)^2)
     * @param quantumNumber The quantum number n (0 = ground state, 1 = first excited, etc.)
     * @param normalized If true, returns properly normalized wavefunction; if false, returns unnormalized
     */
    public HarmonicOscillatorExactSolution(double r0, double alpha, int quantumNumber, boolean normalized) {
        this.r0 = r0;
        this.alpha = alpha;
        this.quantumNumber = quantumNumber;
        this.normalized = normalized;
    }
    
    /**
     * Returns the quantum number (n) of this state.
     * 
     * @return the quantum number
     */
    public int getQuantumNumber() {
        return quantumNumber;
    }
    
    /**
     * Returns the energy of the n-th state.
     * For harmonic oscillator: E_n = (n + 1/2) in the units used (mass = 2, omega = 1)
     * Actually with our potential V = alpha * (r - r0)^2 and mass m = 2:
     * omega = sqrt(2*alpha/m) = sqrt(alpha)
     * E_n = (n + 1/2) * omega = (n + 1/2) * sqrt(alpha)
     * 
     * @return the energy eigenvalue
     */
    public double getEnergy() {
        double omega = Math.sqrt(alpha);
        return (quantumNumber + 0.5) * omega;
    }
    
    /**
     * Evaluates the exact wavefunction at position r.
     * For backward compatibility with ground state (n=0), uses the original formula exp(-alpha*(r-r0)^2).
     * For excited states (n>0), uses the proper harmonic oscillator form: H_n(x) * exp(-x^2/2)
     * 
     * @param r the position
     * @return the wavefunction value ψ(r)
     */
    public double evaluate(double r) {
        // For ground state (n=0), use the original backward-compatible formula
        // Original: exp(-(r - 10)^2) with alpha=1, r0=10
        if (quantumNumber == 0) {
            return Math.exp(-alpha * (r - r0) * (r - r0));
        }
        
        // For excited states (n>0), use proper harmonic oscillator form
        // Reduced coordinate
        double x = Math.sqrt(alpha) * (r - r0);
        
        // Hermite polynomial
        double hermite = hermiteH(quantumNumber, x);
        
        // Gaussian factor: exp(-x^2/2)
        double gaussian = Math.exp(-x * x / 2.0);
        
        // Apply normalization if requested
        if (normalized) {
            double normalization = normalizationConstant(quantumNumber);
            return normalization * hermite * gaussian;
        } else {
            // Return unnormalized (just the Hermite * Gaussian)
            return hermite * gaussian;
        }
    }
    
    /**
     * Evaluates the exact wavefunction at position r using the proper harmonic oscillator form.
     * This method always uses exp(-x^2/2) regardless of quantum number.
     * 
     * @param r the position
     * @return the wavefunction value ψ(r)
     */
    public double evaluateProper(double r) {
        // Reduced coordinate
        double x = Math.sqrt(alpha) * (r - r0);
        
        // Hermite polynomial
        double hermite = hermiteH(quantumNumber, x);
        
        // Gaussian factor (proper form)
        double gaussian = Math.exp(-x * x / 2);
        
        // Apply normalization if requested
        if (normalized) {
            double normalization = normalizationConstant(quantumNumber);
            return normalization * hermite * gaussian;
        } else {
            // Return unnormalized (just the Hermite * Gaussian)
            return hermite * gaussian;
        }
    }
    
    /**
     * Evaluates the exact wavefunction at position r with normalization.
     * 
     * @param r the position
     * @return the normalized wavefunction value ψ(r)
     */
    public double evaluateNormalized(double r) {
        // Reduced coordinate
        double x = Math.sqrt(alpha) * (r - r0);
        
        // Hermite polynomial
        double hermite = hermiteH(quantumNumber, x);
        
        // Gaussian factor
        double gaussian = Math.exp(-x * x / 2);
        
        // Normalization constant
        double normalization = normalizationConstant(quantumNumber);
        
        return normalization * hermite * gaussian;
    }
    
    /**
     * Evaluates the unnormalized (raw) wavefunction at position r.
     * This returns H_n(x) * exp(-x^2/2) without the normalization constant.
     * 
     * @param r the position
     * @return the unnormalized wavefunction value
     */
    public double evaluateUnnormalized(double r) {
        // Reduced coordinate
        double x = Math.sqrt(alpha) * (r - r0);
        
        // Hermite polynomial
        double hermite = hermiteH(quantumNumber, x);
        
        // Gaussian factor
        double gaussian = Math.exp(-x * x / 2);
        
        return hermite * gaussian;
    }
    
    /**
     * Computes the factorial n! using iterative multiplication.
     * 
     * @param n the non-negative integer
     * @return n!
     */
    private double factorial(int n) {
        double result = 1.0;
        for (int i = 2; i <= n; i++) {
            result *= i;
        }
        return result;
    }
    
    /**
     * Computes the normalization constant N_n = 1 / sqrt(2^n * n! * sqrt(pi))
     * 
     * @param n the quantum number
     * @return the normalization constant
     */
    private double normalizationConstant(int n) {
        // N_n = 1 / sqrt(2^n * n! * sqrt(pi))
        double factorial = factorial(n);
        double normalizationFactor = Math.sqrt(Math.pow(2, n) * factorial * Math.sqrt(Math.PI));
        return 1.0 / normalizationFactor;
    }
    
    /**
     * Computes the Hermite polynomial H_n(x) using recurrence relation:
     * H_0(x) = 1
     * H_1(x) = 2x
     * H_{n+1}(x) = 2x*H_n(x) - 2n*H_{n-1}(x)
     * 
     * @param n the order of the Hermite polynomial
     * @param x the argument
     * @return H_n(x)
     */
    private double hermiteH(int n, double x) {
        if (n == 0) {
            return 1.0;
        } else if (n == 1) {
            return 2.0 * x;
        } else {
            // Use recurrence relation
            double h0 = 1.0;
            double h1 = 2.0 * x;
            double h2 = 0.0;
            
            for (int i = 1; i < n; i++) {
                h2 = 2.0 * x * h1 - 2.0 * i * h0;
                h0 = h1;
                h1 = h2;
            }
            
            return h1;
        }
    }
    
    /**
     * Factory method to create ground state (n=0)
     * 
     * @param r0 the center of the potential
     * @param alpha the spring constant
     * @return the ground state exact solution
     */
    public static HarmonicOscillatorExactSolution groundState(double r0, double alpha) {
        return new HarmonicOscillatorExactSolution(r0, alpha, 0);
    }
    
    /**
     * Factory method to create the n-th excited state
     * 
     * @param r0 the center of the potential
     * @param alpha the spring constant
     * @param n the quantum number (0 = ground state, 1 = first excited, etc.)
     * @return the n-th excited state exact solution
     */
    public static HarmonicOscillatorExactSolution excitedState(double r0, double alpha, int n) {
        return new HarmonicOscillatorExactSolution(r0, alpha, n);
    }
}
