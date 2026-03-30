package schrodinger.integrator.pt_correction;

import schrodinger.QuantumLevel;

/**
 * Interface for computing perturbative corrections to the energy.
 * Each correction implementation is specific to a particular integration method.
 */
public interface PTCorrector {
    
    /**
     * Computes and stores the perturbative correction to the energy in the given quantum state.
     * 
     * @param level The quantum state for which to compute the correction
     */
    void computeAndSet(QuantumLevel level);
}