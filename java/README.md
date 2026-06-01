# Schrödinger1D — Java Solver

## Project Summary

Numerical solver for the **1D time-independent Schrödinger equation** (TISE). Given a potential V(r), a particle mass, and a computational grid, it finds eigenvalues (bound-state energies) and eigenfunctions (wavefunctions) using the **shooting method**: propagate ψ from both classical turning points toward a matching interior point, then drive the log-derivative mismatch to zero via bisection or a bisection+regula-falsi hybrid.

The codebase is organized around four orthogonal concerns:
- **Grid**: coordinate mappings from physical r to a uniform computational y
- **Potential**: evaluation of V(r)
- **Integrator**: numerical propagation of ψ (17 interchangeable implementations)
- **Solver**: eigenvalue search strategy

Language: pure Java, no external libraries. Build: IntelliJ IDEA (`.iml` file at workspace root). Tests live in `java/test/` (unit) and `java/test-integration/` (integration).

---

## The Equation Being Solved

The 1D TISE in physical coordinates:

$$-\frac{\hbar^2}{2m}\frac{d^2\psi}{dr^2} + V(r)\psi = E\psi$$

After substituting the coordinate mapping $r \leftrightarrow y$ and defining

$$\tilde{Q}(y) = g(y)^2\bigl[Q(r(y))\bigr] - F(y)$$

$$Q(r) = \frac{2m}{\hbar^2}(E - V(r))$$

the equation becomes the standard Bessel/oscillator normal form

$$\psi''(y) = -\tilde{Q}(y)\,\psi(y)$$

which every integrator solves. Here $g(y) = dr/dy$ is the Jacobian of the mapping and $F(y)$ is a correction term that depends on the mapping strategy. With a uniform grid in y, step size h is constant, which simplifies all recurrences.

All energies are in **Hartree**, distances in **Bohr** (unless noted).

---

## Package Tree

```
java/src/
├── ThreePointProcedural.java          # Standalone pedagogical demo (not integrated)
└── schrodinger/
    ├── Main.java                      # Entry point / usage example
    ├── PhysicalConstants.java         # CODATA 2022 unit conversions
    ├── SchrodingerSystem.java         # Central system object (potential + grid + mass)
    ├── QuantumLevel.java              # Single quantum state (energy + wavefunction)
    ├── FractionalGridCache.java       # Tabulated Q-tilde at fractional indices
    ├── OutputManager.java             # Thread-safe file-based logging
    ├── grid/
    │   ├── Grid.java                  # Computational grid with coordinate transforms
    │   ├── GridFactory.java           # Factory for uniform/log/sqrt/Surkus grids
    │   ├── MappingStrategy.java       # Interface: r ↔ y, g(y), F(y)
    │   ├── MappingLogarithmic.java    # y = r_ref·ln(1 + r/r_ref)
    │   ├── MappingSqrt.java           # y = 2r_ref·(√(1+r/r_ref) − 1)
    │   └── MappingSurkus.java         # y = [(r/r_ref)^α − 1] / [(r/r_ref)^α + 1]
    ├── potential/
    │   ├── PhysicalPotential.java     # Interface: value(r) → Hartree
    │   ├── PhysicalPotentialHarmonic.java
    │   ├── PhysicalPotentialMorse.java
    │   └── PhysicalPotentialLennardJones.java
    ├── integrator/
    │   ├── Integrator.java            # Interface for all propagators
    │   ├── IntegratorFactory.java     # Factory for all 17 integrators
    │   ├── multi_step/
    │   │   ├── Verlet.java            # 2nd order, 2-point history
    │   │   ├── Cowell5.java           # 5th order implicit, 4-point history
    │   │   ├── Cowell6.java           # 6th order implicit, 5-point history
    │   │   ├── Cowell8.java           # 8th order implicit, 7-point history
    │   │   ├── Stormer8.java          # 8th order explicit, 8-point history
    │   │   ├── Obrechkoff6.java       # 6th order, uses Q' and Q''
    │   │   └── numerovlike/
    │   │       ├── Numerov.java               # 4th order classic
    │   │       ├── ExponentiallyFittedAbstract.java  # Base for EF methods
    │   │       ├── EFN.java                   # Exponentially fitted Numerov
    │   │       └── EFNFixedBeta.java          # EF Numerov, β=1/12
    │   ├── predcorr/
    │   │   ├── PCBase.java            # P(EC)^N base framework
    │   │   ├── PC6.java               # 6th order, predictor=EFN
    │   │   ├── PC8Abstract.java       # 8th order corrector base
    │   │   ├── PC8i1.java             # PC8, 1 correction iteration
    │   │   └── PC8i2.java             # PC8, 2 correction iterations
    │   ├── one_step/
    │   │   ├── RK45DP.java            # 5th order Dormand-Prince
    │   │   ├── RKN4.java              # 4th order Runge-Kutta-Nyström
    │   │   ├── CFMagnusAbstract.java  # Commutator-free Magnus base
    │   │   ├── CFMagnus4.java         # 4th order Magnus, 2 exponentials
    │   │   ├── CFMagnus6e4.java       # 6th order Magnus, 4 exponentials
    │   │   ├── CFMagnus6e5Opt.java    # 6th order Magnus, 5 exponentials (optimized)
    │   │   └── CFMagnus8.java         # 8th order Magnus, 11 exponentials (default best)
    │   └── pt_correction/
    │       ├── PTCorrector.java       # Interface: post-hoc energy correction
    │       ├── NumerovPTCorrector.java
    │       └── VerletPTCorrector.java
    └── solver/
        ├── ShootingSolver.java        # Main eigenvalue finder
        ├── EigenvalueBounds.java      # Energy bracket manager per quantum number
        └── RefinementStrategy.java    # Enum: BISECTION_ONLY | BISECTION_THEN_BIDIRECTIONAL
```

---

## Core Classes

### `Main` (`schrodinger.Main`)

Entry point demonstrating the canonical usage pattern. Creates a harmonic oscillator, builds a uniform grid, constructs a `SchrodingerSystem`, picks an integrator, and calls `ShootingSolver.findEigenvaluesUpTo()`. Contains only formatting helpers (`padFloat`, `padInt`, `padding`) and `main()`. Intended as a living usage example.

---

### `PhysicalConstants` (`schrodinger.PhysicalConstants`)

Static constants and conversion methods based on CODATA 2022. No instances.

| Constant | Value | Meaning |
|---|---|---|
| `HARTREE_TO_INVERSE_CM` | 219 474.631 363 14 | Hartree → cm⁻¹ |
| `BOHR_TO_ANG` | 0.529 177 210 544 | Bohr → Ångström |
| `UMA_TO_ELECTRON_MASS` | 1 822.888 486 278 | u → mₑ |

Key methods: `toInverseCm(double hartree)`, `toHartree(double cm1)`.

---

### `SchrodingerSystem` (`schrodinger.SchrodingerSystem`)

The central object passed through the entire solver pipeline. Holds:

| Field | Type | Purpose |
|---|---|---|
| `physicalPotential` | `PhysicalPotential` | V(r) |
| `mass` | `double` | Particle mass (atomic units) |
| `grid` | `Grid` | Coordinate grid |
| `qMin` | `double` | Cap floor for Q-tilde (prevents overflow in classically forbidden regions) |
| `UTildeMinimumGridIndex` | `double` | Grid index of potential minimum |
| `UTildeMinimumGridValue` | `double` | Value of effective potential at its minimum |
| `uMaxRight` | `double` | Max effective potential to the right of minimum |
| `maxStepSizeAllowedRegion` | `double` | Largest h still inside the classically allowed region |
| `cacheUTilde` | `FractionalGridCache` | Optional: precomputed Q-tilde at fractional indices |

Key methods:

- `U(double r)` — raw V(r) in Hartree
- `UTildeAtGridPoint(double i)` — effective potential $\tilde{U}(y_i)$ including mapping correction; uses cache if initialized
- `Ucapped()` — returns `qMin`-capped value (used when `isCapPotential()` is true)
- `estimateEnergyScale()` — typical energy scale derived from grid spacing and Q-tilde; used for initial bracket estimation
- `findLeftmostInversionGridpointY(double energy)` / `findRightmostInversionGridpointY(double energy)` — binary search for classical turning points
- `initializeCache(double[] offsets)` — builds `FractionalGridCache` for the supplied fractional offsets
- `isCachingUTilde()`, `setCachingUTilde(boolean)` — toggle caching

---

### `QuantumLevel` (`schrodinger.QuantumLevel`)

Represents one quantum state throughout the solver lifecycle. Passed into integrators for propagation and returned to the caller after convergence.

| Field | Type | Purpose |
|---|---|---|
| `system` | `SchrodingerSystem` | Back-reference to system |
| `energy` | `double` | Current trial or converged eigenvalue |
| `lowerBound`, `upperBound` | `double` | Energy bracket from bisection |
| `nodesLower`, `nodesUpper` | `int` | Node counts at bracket endpoints |
| `psi` | `double[]` | Wavefunction values at grid points |
| `currentPsiPrime` | `double[]` | Wavefunction derivative (needed by one-step methods) |
| `perturbativeCorrectionToEnergy` | `double` | Post-hoc PT energy correction |
| `convergenceInfo` | `List<ConvergenceInfo>` | Per-stage iteration counts for diagnostics |
| `isCapPotential` | `boolean` | Whether Q-tilde capping is active |

Key methods:

- `Q(double r)` — $Q(r) = 2m(E - V(r))/\hbar^2$
- `QTilde(double y)` / `QTildeAtGridPoint(double i)` — transformed Q (with capping if `isCapPotential`)
- `QTildePrimeAtGridPoint(double i)` / `QTildeDoublePrimeAtGridPoint(double i)` — numerical derivatives (used by Obrechkoff6)
- `normalizePsi()` — normalizes $\psi$ to unit L² norm; returns the norm before normalization
- `countNodes()` — counts zero crossings in `psi`
- `maximumStepSize()` — returns the largest h consistent with accuracy for the current energy
- `verifyStepSize()` — logs a warning if the grid step exceeds `maximumStepSize()`
- `countTotalScans()` — sums `convergenceInfo` iteration counts

Inner class `ConvergenceInfo(String stage, int iterations)` records how many propagations each solver stage consumed.

---

### `FractionalGridCache` (`schrodinger.FractionalGridCache`)

Precomputes and stores a function's values at integer and fractional grid indices (e.g., n, n+0.5, n+0.25) so integrators can retrieve them in O(1) instead of recomputing per step.

Construction: `new FractionalGridCache(int numGridPoints, double[] offsets, DoubleUnaryOperator function)`.

- `get(double i)` — returns cached value at (possibly fractional) index i; resolves which offset slot to use automatically
- `getByFlatIndex(int flatIndex)` — direct access into the flattened backing array
- `computeCacheVectorIndex(double i)` — maps a fractional index to its flat array position
- Exposes hit/miss counters (`nOfCacheHits`, `nOfCacheMisses`) for diagnostics

Integrators declare their required offsets via `Integrator.getFractionalOffsets()`; `SchrodingerSystem.initializeCache()` builds the cache for those offsets.

---

### `OutputManager` (`schrodinger.OutputManager`)

Thread-safe, timestamped file logging. All methods are static and synchronized.

- `initCommonOutputFile(String filePath)` — opens the shared log file
- `write(String msg)` / `writeBlankLine()` / `writeData(String line)` — write to shared file
- `createLevelWriter(String filePath)` → `LevelWriter` — per-eigenstate log file

Inner class `LevelWriter` wraps a per-level `PrintWriter`: `log(msg)` (timestamped), `raw(line)` (verbatim), `close()`.

---

## Grid Package (`schrodinger.grid`)

### `Grid`

Encapsulates the uniform computational grid in y, with bidirectional transforms to physical r.

| Field | Type | Purpose |
|---|---|---|
| `rMin`, `rMax` | `double` | Physical coordinate bounds |
| `yMin`, `yMax` | `double` | Transformed coordinate bounds |
| `stepSizeYCoordinate` | `double` | Uniform step h in y |
| `numberOfPoints` | `int` | Grid size N |
| `mappingStrategy` | `MappingStrategy` | The active r ↔ y transform |

Key methods:

- `yAtGridPoint(double i)` — `yMin + i·h`
- `rAtGridPoint(double i)` — `r(yAtGridPoint(i))`
- `y(double r)` / `r(double y)` — delegate to `mappingStrategy`
- `g(double y)` / `gAtGridPoint(double i)` — `dr/dy` Jacobian
- `F(double y)` — mapping correction term for Q-tilde
- `getFirstYValue()`, `getLastYValue()`, `getStepSizeYCoordinate()`, `getNumberOfPoints()`, `getMappingStrategy()`

---

### `GridFactory`

Static factory; all methods return a new `Grid`.

| Method | Mapping |
|---|---|
| `generateUniformGrid(rMin, rMax, N)` | Identity (no transform) |
| `generateLogarithmicGrid(rMin, rMax, N, rRef)` | Logarithmic |
| `generateSqrtGrid(rMin, rMax, N, rRef)` | Square root |
| `generateSurkusGrid(rMin, rMax, N, rRef, alpha)` | Surkus |

---

### `MappingStrategy` (interface)

Default implementations are the identity map. Implementations override some or all of:

- `y(double r)` — physical → computational
- `r(double y)` — computational → physical
- `g(double y)` — `dr/dy`
- `F(double y)` — $g''/(2g) - (g')^2/(4g^2)$ correction to Q-tilde

**`MappingLogarithmic`**: $y = r_\text{ref}\ln(1 + r/r_\text{ref})$, $g = e^{y/r_\text{ref}}$.
Useful for long-range potentials; compresses large r while expanding small r.

**`MappingSqrt`**: $y = 2r_\text{ref}(\sqrt{1 + r/r_\text{ref}} - 1)$, $g = (1 + y/(2r_\text{ref}))$.
Smoother compression than logarithmic.

**`MappingSurkus`**: $y = \frac{(r/r_\text{ref})^\alpha - 1}{(r/r_\text{ref})^\alpha + 1}$, range of y is $(-1, 1)$.
Flexible power-law mapping; α controls the rate of compression. Useful for potentials with both short-range repulsion and long-range tails (e.g., Lennard-Jones).

---

## Potential Package (`schrodinger.potential`)

### `PhysicalPotential` (interface)

Single method: `double value(double r)` returning V(r) in Hartree.

### `PhysicalPotentialHarmonic`

$$V(r) = \alpha\,(r - r_0)^2$$

Fields: `r0` (equilibrium), `alpha` (force constant). Used for testing and demos.

### `PhysicalPotentialMorse`

$$V(r) = D_e\bigl[1 - e^{-a(r - r_\text{min})}\bigr]^2$$

Fields: `rmin`, `a` (range parameter), `De` (well depth). Models realistic diatomic bonds with correct dissociation behavior.

### `PhysicalPotentialLennardJones`

$$V(r) = 4\varepsilon\left[\left(\frac{\sigma}{r}\right)^n - \left(\frac{\sigma}{r}\right)^{n/2}\right] + \varepsilon$$

Fields: `wellDepth` (ε), `sigma`, `exponent` (n, default 6). Models van der Waals / noble-gas interactions. The `+ε` shift sets the dissociation limit to zero.

---

## Integrator Package (`schrodinger.integrator`)

### `Integrator` (interface)

All propagators implement:

```java
double propagate(
    double[] psi,                    // wavefunction history buffer (length ≥ minHistoryLength())
    double[] currentPsiPrime,        // derivative buffer (for one-step methods; else null)
    int n,                           // current grid index
    double step,                     // h
    DoubleUnaryOperator qTilde,      // Q̃(i) as function of fractional index
    DoubleUnaryOperator qTildePrime, // dQ̃/di
    DoubleUnaryOperator qTildeDoublePrime, // d²Q̃/di²
    Direction direction              // FORWARD or BACKWARD
);
```

Returns the new value of ψ at n+1 (or n−1 for `BACKWARD`). The caller is responsible for maintaining the `psi` ring buffer.

Additional contract methods:

| Method | Meaning |
|---|---|
| `minHistoryLength()` | Minimum entries in `psi` buffer |
| `globalConvergenceOrder()` | Formal order p (error ~ hᵖ) |
| `needsPotentialCapping()` | Must Q-tilde be capped? (multi-step methods: yes) |
| `getFractionalOffsets()` | Fractional indices where Q-tilde must be pre-tabulated |
| `getPertubativeCorrector()` | Associated `PTCorrector`, or null |

`Direction` enum: `FORWARD` (+1), `BACKWARD` (−1).

---

### `IntegratorFactory`

Static factory returning every integrator by name. Notable:

- `getBestOneStepIntegrator()` → `CFMagnus8` (used internally by `ShootingSolver` for verification passes)
- `getAll()` → `List<Integrator>` of all 17 integrators (useful for benchmark loops)

---

### Integrator Overview Table

| Class | Order | History | Capping | Category |
|---|---|---|---|---|
| `Verlet` | 2 | 2 | No | multi-step |
| `Numerov` | 4 | 2 | Yes | Numerov-like |
| `EFNFixedBeta` | 4 | 2 | Yes | Numerov-like |
| `EFN` | 4 | 2 | No | Numerov-like |
| `Cowell5` | 5 | 4 | Yes | multi-step |
| `PC6` | 6 | 3 | Yes | predictor-corrector |
| `Cowell6` | 6 | 5 | Yes | multi-step |
| `Obrechkoff6` | 6 | 2 | Yes | multi-step |
| `RKN4` | 4 | 1 | No | one-step |
| `CFMagnus4` | 4 | 1 | No | one-step |
| `RK45DP` | 5 | 1 | No | one-step |
| `CFMagnus6e4` | 6 | 1 | No | one-step |
| `CFMagnus6e5Opt` | 6 | 1 | No | one-step |
| `Cowell8` | 8 | 7 | Yes | multi-step |
| `Stormer8` | 8 | 8 | No | multi-step |
| `PC8i1` | 8 | 4 | Yes | predictor-corrector |
| `PC8i2` | 8 | 4 | Yes | predictor-corrector |
| `CFMagnus8` | 8 | 1 | No | one-step |

---

### Multi-Step Methods (`multi_step/`)

All store a ring buffer of past ψ values. "Capping required" means the method diverges if Q-tilde is allowed to grow unboundedly in the forbidden region.

**`Verlet`** — $\psi_{n+1} = 2\psi_n - \psi_{n-1} - h^2\tilde{Q}_n\psi_n$. Simplest; has a `VerletPTCorrector`.

**`Cowell5/6/8`** — implicit linear multi-step (Cowell–Numerov family). The implicit value $\psi_{n+1}$ appears on both sides but can be solved in closed form because the equation is linear. Coefficients are fixed rational numbers. Order 5/6/8 with 4/5/7 history points respectively.

**`Stormer8`** — explicit 8th-order Störmer formula; 8-point history of Q-tilde values. No implicit step.

**`Obrechkoff6`** — 6th order; uses $\tilde{Q}'$ and $\tilde{Q}''$ in addition to $\tilde{Q}$. Only 2 history points needed despite higher order. Requires `QTildePrimeAtGridPoint` and `QTildeDoublePrimeAtGridPoint`.

---

### Numerov-like Methods (`multi_step/numerovlike/`)

**`Numerov`** — classic 4th-order formula:

$$\psi_{n+1} = \frac{2\psi_n(1 - \tfrac{5}{6}h^2\tilde{Q}_n) - \psi_{n-1}(1 + \tfrac{1}{12}h^2\tilde{Q}_{n-1})}{1 + \tfrac{1}{12}h^2\tilde{Q}_{n+1}}$$

Has a `NumerovPTCorrector`. The reference method for benchmarking.

**`ExponentiallyFittedAbstract`** — base class for methods that adapt their coefficients to the local oscillatory vs. exponential character of the solution (determined by the sign and magnitude of $\tilde{Q}$). Subclasses implement `getBeta(Z)` and `getGamma(Z, beta)` where $Z = h^2\tilde{Q}$.

**`EFN`** (Exponentially Fitted Numerov) — exact fitting: coefficients derived from $\sin(\sqrt{Z}/2)$ (oscillatory) or $\sinh(\sqrt{|Z|}/2)$ (exponential) with Taylor fallback for $|Z| < 0.01$.

**`EFNFixedBeta`** — simplified variant: $\beta = 1/12$ (constant), $\gamma$ still adapted. Slightly more robust with capping.

---

### Predictor-Corrector Methods (`predcorr/`)

Implement P(EC)^N schemes: predict ψ_{n+1} with a lower-order integrator, then correct N times with a higher-order symmetric stencil.

**`PCBase`** — base framework. Holds a `predictor: Integrator` injected at construction. Provides `predictAhead()` which builds the prediction chain over multiple steps (needed to fill the corrector's symmetric stencil).

**`PC6`** — 6th order. Predictor: `EFN` (4th order). Corrector: 5-point symmetric stencil {−2,−1,0,+1,+2}:

$$\psi_n = \sum_{k=-2}^{2} \beta_k h^2 \tilde{Q}_{n+k}\psi_{n+k} + \psi_{n-2} + \psi_{n+2} - ...$$

**`PC8Abstract`** — 8th order corrector base: 7-point stencil {−3,...,+3}. Subclasses differ in iteration count.

**`PC8i1`** — 1 correction iteration after initial prediction.

**`PC8i2`** — 2 correction iterations; more accurate but slower.

---

### One-Step Runge-Kutta Methods (`one_step/`)

Require history length 1 (only current point). Propagate both ψ and ψ' simultaneously. Never require capping. Well-suited for rapidly varying potentials.

**`RKN4`** — 4th order Runge-Kutta-Nyström for second-order ODEs. Uses half-step node (fractional offset 0.5). Fewer stages than general RK.

**`RK45DP`** — 5th order Dormand-Prince. Standard 7-stage tableau. 7 function evaluations per step. Fractional offsets: {0, 1−c5, c2, c3, 1−c3, c4, c5} (Dormand-Prince nodes).

---

### Commutator-Free Magnus Methods (`one_step/`)

All extend `CFMagnusAbstract`. These treat the ODE as a linear system $\psi'' = -\tilde{Q}\psi$ and apply the Magnus exponential series approximated at Gauss-Legendre quadrature nodes. The result is a product of matrix exponentials (2×2 rotation/hyperbolic matrices) that automatically adapts between oscillatory and exponential regimes.

**`CFMagnusAbstract`** — base. Precomputes weight matrices from Legendre polynomial collocation. The `propagate` method evaluates $\tilde{Q}$ at the Gauss-Legendre nodes within the step, constructs the exponentials (cos/cosh, sin/sinh depending on sign of $\tilde{Q}$), and applies them in reverse order.

Static utility `computeW(...)` builds the weight matrices at class initialization from the method's F-matrix and Gauss-Legendre data.

| Class | Order | GL nodes | Exponentials |
|---|---|---|---|
| `CFMagnus4` | 4 | 2 | 2 |
| `CFMagnus6e4` | 6 | 3 | 4 |
| `CFMagnus6e5Opt` | 6 | 3 | 5 (optimized stability) |
| `CFMagnus8` | 8 | 4 | 11 |

`CFMagnus8` is the default "best" integrator (`IntegratorFactory.getBestOneStepIntegrator()`). It is significantly slower per step (~1800 ns/point) than multi-step methods (~20–60 ns/point) but requires no bootstrapping and handles non-smooth potentials well.

---

## Perturbative Correctors (`integrator/pt_correction/`)

Applied after convergence to improve the energy estimate without additional propagations.

**`PTCorrector`** (interface) — `computeAndSet(QuantumLevel level)`: reads `level.psi` and modifies `level.perturbativeCorrectionToEnergy`.

**`NumerovPTCorrector`** — correction from finite-difference second derivatives of Q-tilde:

$$\Delta E = \frac{h^3}{480\,m}\sum_i \bigl(\tilde{Q}'_{i+1}\psi_{i+1} - \tilde{Q}'_{i-1}\psi_{i-1}\bigr)^2$$

**`VerletPTCorrector`** — simpler point-value correction:

$$\Delta E = \frac{h^3}{24\,m}\sum_i \bigl(\tilde{Q}_i\,\psi_i\bigr)^2$$

---

## Solver Package (`schrodinger.solver`)

### `ShootingSolver`

The main driver. Constructor: `ShootingSolver(SchrodingerSystem system, Integrator integrator)`.

Key constants:
- `TARGET_ABSOLUTE_ERROR = 1e-13` (Hartree) — bisection stopping criterion
- `MAXIMUM_NUMBER_OF_BISECTIONS = 60` — hard limit per eigenvalue
- `PSI_MAX = 1e140` — threshold for mid-propagation rescaling (prevents overflow in forbidden regions)

#### Algorithm

1. **`findInitialEnergyBracket(int nNodes)`** — scans energies to find an interval $[E_\text{low}, E_\text{high}]$ where the propagated node count brackets the target. Uses `EigenvalueBounds` to reuse information from previously found states.

2. **`findEigenvalueByBisection(QuantumLevel, int matchingIndex)`** — pure bisection: trial energy = midpoint, propagate from both ends to `matchingIndex`, count nodes, update bracket until `upperBound − lowerBound < TARGET_ABSOLUTE_ERROR`.

3. **`findEigenvalueByHybridMethod(QuantumLevel, int matchingIndex)`** — same bisection loop but switches to regula falsi (using derivative mismatch as the function) once the bracket is tight enough. Converges faster near the root.

4. **`computeDerivativeMismatch(QuantumLevel, int matchingIndex)`** — propagates ψ from left and right, computes log-derivative $(ψ'/ψ)$ from each side at `matchingIndex`, returns their difference. This is the residual driven to zero.

5. **`findMatchingIndex(double energy)`** — locates the classical turning point (or a fixed interior point) as the matching location.

6. After convergence: calls `PTCorrector.computeAndSet(level)` if the integrator provides one, then `level.normalizePsi()`.

#### Public API

- `findEigenvaluesUpTo(int maxN, RefinementStrategy)` → `List<QuantumLevel>` — find states v=0…maxN
- `findEigenvalue(int nNodes)` → `QuantumLevel` — find single state with default strategy
- `findEigenvalue(int nNodes, RefinementStrategy)` → `QuantumLevel` — find single state

---

### `EigenvalueBounds`

Tracks per-quantum-number energy brackets across successive `findEigenvalue` calls so that the bracket search for state v+1 can reuse the upper bound of state v.

- `updateBounds(double energy, int nNodes)` — updates `lowerBounds[nNodes]` and `upperBounds[nNodes]`
- `getLowerBound(int v)` / `getUpperBound(int v)`
- `isLowerBoundDefined(int v)` / `isUpperBoundDefined(int v)`
- `getLowerNodes(int v)` / `getUpperNodes(int v)`

---

### `RefinementStrategy` (enum)

| Value | Behaviour |
|---|---|
| `BISECTION_ONLY` | Pure bisection; always converges, O(log(1/ε)) iterations |
| `BISECTION_THEN_BIDIRECTIONAL` | Bisection until bracket is narrow, then bidirectional derivative-matching (faster near root) |

---

## Data Flow

```
PhysicalPotential  ─┐
Grid (GridFactory) ─┤
mass               ─┴─▶  SchrodingerSystem
                              │
                              ▼
                    Integrator (IntegratorFactory)
                              │
                              ▼
                       ShootingSolver
                              │
                   ┌──────────┴──────────┐
                   │  findInitialBracket  │
                   │  (EigenvalueBounds)  │
                   └──────────┬──────────┘
                              │
                   ┌──────────┴──────────┐
                   │   bisection loop    │
                   │  propagate forward  │
                   │  propagate backward │
                   │  count nodes / mismatch │
                   └──────────┬──────────┘
                              │
                   ┌──────────┴──────────┐
                   │   PTCorrector       │  (if integrator has one)
                   │   normalizePsi()    │
                   └──────────┬──────────┘
                              │
                              ▼
                         QuantumLevel
                  (energy, psi[], convergenceInfo)
```

Propagation inside the bisection loop calls `Integrator.propagate()` at each grid step, which reads Q-tilde via the `DoubleUnaryOperator` passed in (backed by `FractionalGridCache` if caching is active).

---

## Design Patterns

| Pattern | Where |
|---|---|
| **Strategy** | `Integrator` (17 impls), `MappingStrategy` (4 impls), `PhysicalPotential` (4 impls) |
| **Factory** | `IntegratorFactory`, `GridFactory` |
| **Template Method** | `ExponentiallyFittedAbstract` (subclasses supply β, γ), `CFMagnusAbstract` (subclasses supply node/weight tables), `PC8Abstract` (subclasses supply iteration count) |
| **Caching** | `FractionalGridCache` precomputes Q-tilde at offsets declared by integrators |
| **Builder / Fluent** | `SchrodingerSystem` constructor parameters |

---

## Build & Tests

**Build system**: IntelliJ IDEA only (`Schodinger1D.iml` at workspace root). No Maven or Gradle. No external dependencies — pure Java standard library.

**Unit tests** — `java/test/schrodinger/`:
- `FractionalGridCacheTest`, `QuantumLevelTest`
- `grid/`, `integrator/` subdirectories with per-class tests

**Integration tests** — `java/test-integration/schrodinger/`:
- `IntegrationTest`, `ConvergenceTest`, `DissociationTest`
- `morse/`, `lennard_jones/` subdirectory tests
- `Utils` — shared test utilities
- Reference data in `java/test/resources/`: `AiryAi_reference_values_3025.csv`, `MathieuS_reference_values_3025.csv`

---

## Fortran Companion Files (Reference / Benchmark)

The workspace root contains Fortran code that is **not part of the Java solver** but serves as reference material:

| Path | Purpose |
|---|---|
| `benchS.FOR`, `benchD.FOR` | Single- and double-precision benchmarks of Numerov-type methods |
| `benchmark/` | Additional Fortran benchmark variants (`FINAL.FOR`, `for_bench.f90`, etc.) |
| `level16/` | LEVEL 16 Fortran program (Le Roy) — established Schrödinger solver used for validation |
| `roundoff/` | Fortran round-off error studies comparing integration methods |
| `oofortran/` | Object-oriented Fortran prototype |
| `ThreePoint_*.f90` | Fortran three-point propagation reference implementations |

The `level16/` directory includes the full LEVEL 16 documentation (`level16documentation_01_index.txt` through `_05`) and sample I/O files used to validate the Java solver's eigenvalues.

---

## Standalone Demo

**`java/src/ThreePointProcedural.java`** — a self-contained, dependency-free demonstration of three-point propagation on a harmonic oscillator $V(x) = x^2$. Uses hardcoded bisection and tests 20 step sizes from h=0.1 down to h=10⁻⁶·²⁵. Not integrated with the solver framework; read it to understand the bare algorithm before studying `ShootingSolver`.

---

## Known Limitations / Roadmap (from `java/docs/TODO.md`)

- No rotational term (centrifugal barrier) yet
- No input file format (potential and grid parameters are hardcoded in `Main`)
- Missing grid mappings: linear-switchover, arctan, rational
- No polynomial/spline potential type
- Matrix element calculations (Einstein coefficients) not yet implemented
- No Richardson extrapolation post-processor
- Non-adiabatic corrections not yet planned
- Fortran port of selected integrators under consideration
