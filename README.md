# evolving-sheaf-c

**Spectral Gap Dynamics in Evolving Cellular Sheaves**

> When a theorem fails, the failure is data. The spectral gap of an evolving sheaf CAN decrease. This library studies WHEN, HOW FAST, and WHERE it happens.

---

## The Failure That Taught Us Something

Classical sheaf cohomology over graphs gives us the Hodge Laplacian and its spectral gap λ₁ — the smallest nonzero eigenvalue, a measure of how tightly the sheaf glues local data into global consistency. The theorem is clean and reassuring: for a **static** sheaf, the spectral gap is a fixed invariant. It doesn't budge. Time passes, energy flows, but λ₁ sits there like a rock.

Then we let the restriction maps breathe.

When restriction maps become functions of flow energy — when they *evolve* — the spectral gap stops being a constant and starts being a signal. On cycle-4 with linear evolving sheaves and sinusoidal flow, the gap drops from 5.36 to 4.84. That's a 9.7% erosion. Not a rounding error. Not a numerical artifact. A genuine, reproducible decrease in the sheaf's global connectivity.

This wasn't caught by manual testing. It was caught by a **DeepSeek audit** that flagged the decreasing gap as a potential implementation bug. We investigated. It wasn't a bug. The math was correct. The *theorem* was what broke — or rather, the theorem's assumption of static structure no longer held. The audit didn't find a defect. It found a discovery.

**This library is the instrument for studying that discovery.**

## The Story in Brief

The static sheaf theorem says: build your Hodge Laplacian L = D†D from fixed restriction maps, and λ₁ is a topological invariant. Constant. Unchanging. Beautiful.

But real systems aren't static. Distributed consensus protocols have messages in flight. Sensor networks have varying signal strength. Social networks have evolving trust. In all of these, the "restriction maps" — the local-to-local compatibility conditions — depend on time-varying quantities. We model this as **flow energy** E(t) that modulates the restriction maps.

When R(t) = R₀ + α·E(t), something fundamental shifts. The spectral gap becomes a function of time. It oscillates. It erodes. It can recover, then erode again. And crucially, different graph topologies respond *differently*. Expander graphs shrug off the perturbation. Cycle graphs get hammered. Path graphs get destroyed.

This isn't just mathematically interesting. It's *practically* important. If your distributed system's convergence rate depends on the spectral gap, and the gap is shrinking, your system is slowing down and you don't even know it.

## What's Inside

### Three Sheaf Models

The library implements three restriction map regimes:

| Model | Restriction Map | Spectral Gap Behavior |
|-------|----------------|----------------------|
| **Static** | R(t) = R₀ | Constant — the theorem holds ✓ |
| **Linear Evolving** | R(t) = R₀ + α·E(t) | Decreases with energy — the theorem breaks ✗ |
| **Nonlinear Evolving** | R(t) = R₀ · f(E(t)) | Complex dynamics — rich phase structure ✗ |

The nonlinear model supports three activation functions: **sigmoid** (bounded saturation), **tanh** (odd symmetry with sign changes), and **exponential decay** (Gaussian suppression). Each produces qualitatively different gap dynamics.

### Phase Transition Detection

When dλ₁/dt changes sign, something fundamental shifts in the sheaf's connectivity structure. The library tracks these **phase transitions** across trajectories. Static sheaves have zero transitions — as expected. Linear evolving sheaves under sinusoidal flow produce regular, periodic transitions. Pulse flows create sharp, sudden transitions that look almost like first-order phase changes.

The phase transitions aren't just detection artifacts. They mark real moments where the system transitions from "losing coherence" to "regaining coherence" (or vice versa). In a distributed system, these are the moments where convergence rate switches from degrading to improving.

### Stability Maps

The `es_stability_map()` function computes the spectral gap over a grid of (α, time) parameter values. This produces a 2D surface showing exactly where the gap grows, where it shrinks, and where the critical boundaries lie. For a fixed graph topology, you can visualize the entire parameter space of flow-perturbation strength vs. time.

At α = 0, the surface is flat — the static theorem holds. As α increases, valleys form. The topology of those valleys depends on the graph topology.

### Topology Comparison: Expanders Win

The most striking result so far:

| Topology | Relative Gap Change | Interpretation |
|----------|-------------------|----------------|
| **Expander** (3-regular) | ~0.24 | Robust. High Cheeger constant absorbs the perturbation. |
| **Cycle** | ~0.50 | Fragile. Low expansion means the gap gets hammered. |
| **Path** | Worst | Maximally fragile. Single path of failure propagation. |

The message is clear: **if you want your system to be robust against evolving restriction maps, build it like an expander.** High connectivity, good expansion, no bottlenecks. The Cheeger constant isn't just a graph invariant — it's a robustness predictor.

## Build & Run

```bash
# Build and test
make test

# Or manually
gcc -O2 -I include -o test_runner tests/test_evolving_sheaf.c src/evolving_sheaf.c -lm
./test_runner
```

No external dependencies. The eigenvalue solver uses Jacobi rotation — no LAPACK, no BLAS, no nothing. Pure C.

## Quick Example

```c
#include "evolving_sheaf.h"
#include <stdio.h>

int main() {
    // Cycle graph C4
    es_graph *g = es_graph_cycle(4);

    // Linear evolving sheaf: R = 1.0 + 0.5·E(t)
    es_sheaf_config cfg = {
        .model = ES_LINEAR,
        .R0 = 1.0,
        .alpha = 0.5
    };

    // Track gap from t=0 to t=20
    es_gap_trajectory *traj = es_track_gap(
        g, &cfg, es_flow_sinusoidal, NULL,
        0.0, 20.0, 200
    );

    printf("Gap: min=%.3f, max=%.3f, transitions=%d\n",
        traj->min_gap, traj->max_gap, traj->n_transitions);
    printf("Total change: %.3f (relative: %.3f)\n",
        traj->total_change, traj->total_change / traj->max_gap);

    es_gap_trajectory_free(traj);
    es_graph_free(g);
    return 0;
}
```

## Applications

### Self-Healing Distributed Systems

If you can monitor the spectral gap of your system's sheaf Laplacian in real time, you can detect when coherence is degrading. Phase transitions mark the tipping points. A system that detects its gap shrinking can trigger remediation — rerouting, load rebalancing, or topology repair.

### Convergence Rate Prediction

The spectral gap directly controls convergence rates in consensus and distributed optimization algorithms. If the gap is decreasing, your algorithm is slowing down. This library lets you predict *when* that happens and *how bad* it gets for a given perturbation strength.

### Network Robustness Analysis

The expander-vs-cycle comparison generalizes. For any network topology, you can compute how much the spectral gap degrades under flow perturbation. This is a new measure of network robustness that goes beyond static spectral analysis — it measures *dynamic* robustness.

## The Research Frontier

This library is a tool, not a conclusion. The open questions are the point:

**1. Is there a universal bound on gap decrease rate?**

For linear evolving sheaves R(t) = R₀ + α·E(t), does |dλ₁/dt| have a universal upper bound in terms of α and the graph's Cheeger constant? Preliminary evidence suggests O(α · h(G)⁻¹). If true, this would be a genuine *theorem about theorem failure* — a bound on how fast the static theorem can break.

**2. Which topologies resist gap decrease?**

The expander result is suggestive but not proven. Is it the Cheeger constant specifically that matters, or is it something deeper about the graph's spectral geometry? Are there topologies that are *more* robust than expanders? Random graphs with planted expanders?

**3. Is there a gap conservation law?**

When λ₁ decreases, do higher eigenvalues compensate? If the spectral gap erodes, does the "missing" spectral weight appear elsewhere in the spectrum? A conservation law would mean that evolving sheaves don't lose connectivity — they *redistribute* it. This would be a spectral equipartition principle hidden in the Hodge structure.

**4. Phase transitions as critical phenomena**

The sign change in dλ₁/dt marks a phase transition. But what's the order parameter? Is it continuous or first-order? The sigmoid nonlinear model shows particularly rich transition behavior — is this a universality class? Can we classify the phase diagrams of evolving sheaves?

## Architecture

```
include/evolving_sheaf.h    — Public API with full documentation
src/evolving_sheaf.c        — Core: Hodge Laplacian construction + Jacobi eigensolver
tests/test_evolving_sheaf.c — 43 tests: models, transitions, stability, topologies
```

### Eigenvalue Solver

The Hodge Laplacian is symmetric positive semi-definite. We compute all eigenvalues using Jacobi rotation — iterative plane rotations that annihilate off-diagonal entries until the matrix is diagonal. Convergence to ~1e-14 precision. No external dependencies. Handles up to at least n=20 comfortably (tested), and larger with O(n³) patience.

### Flow Energy Models

Four built-in flow generators:
- **Constant**: E(t) = 1.0 — reduces to static
- **Sinusoidal**: E(t) = 1 + 0.5·sin(0.5t + edge·0.3) — smooth periodic perturbation with edge-dependent phase
- **Pulse**: sharp periodic spikes — stress-tests gap response
- **Random walk**: deterministic pseudo-random — tests irregular perturbation patterns

Custom flows via `es_flow_fn` function pointer.

## Test Coverage

43 tests across 10 categories. Every model. Every flow. Edge cases (single edge, large graphs, extreme α). The expander robustness test doesn't just assert — it prints the observed relative changes so you can see the phenomenon yourself.

The most important test might be `test_gap_sum_experiment` — it records trace observations across time points, always passing, because it's not asserting a hypothesis. It's *collecting data* for the gap conservation question.

## The Philosophical Point

Most work on sheaf cohomology assumes static structure. The theorem is elegant: the spectral gap is an invariant, a fixed number that characterizes the sheaf's global connectivity. Clean, simple, publishable.

But nature isn't static. Networks evolve. Trust changes. Signal strength varies. Flows propagate and dissipate. When we let the sheaf model that reality — when restriction maps become functions of energy rather than constants — the spectral gap stops being a number and starts being a *story*. It tells you when the system is losing coherence, when it's recovering, and how much it can take before things fall apart.

The DeepSeek audit flagged a "bug." But the bug was reality. The theorem's assumption was wrong, not the implementation.

**The failure of the theorem IS the research program.**

## License

MIT
