package schrodinger.pt_correction;

import schrodinger.QuantumState;

/**
 * Interface for computing perturbative corrections to the energy.
 * Each correction implementation is specific to a particular integration method.
 */
public interface PTCorrection {
    
    /**
     * Computes and stores the perturbative correction to the energy in the given quantum state.
     * 
     * @param level The quantum state for which to compute the correction
     */
    void compute(QuantumState level);
}