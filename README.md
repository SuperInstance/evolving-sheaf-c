# evolving-sheaf-c

**Spectral Gap Dynamics in Evolving Cellular Sheaves**

> When a theorem *fails*, the failure is more interesting than the theorem.
> The evolving sheaf spectral gap decreases in dynamic cases — and that decrease *is* the signal.

## The Discovery

Classical sheaf cohomology gives us the Hodge Laplacian and its spectral gap λ₁ — a measure of how "connected" a sheaf is over a graph. The theorem says: for a *static* sheaf, the spectral gap is a fixed invariant.

But when restriction maps **evolve with flow energy** — when the sheaf *breathes* — the spectral gap changes. On cycle-4, it drops from 5.36 to 4.84. That's not a bug. That's a **phenomenon**.

This library studies that phenomenon.

## What It Does

### Three Sheaf Models

| Model | Restriction Map | Gap Behavior |
|-------|----------------|--------------|
| **Static** | R(t) = R₀ | Constant (theorem holds ✓) |
| **Linear evolving** | R(t) = R₀ + α·E(t) | Decreases with energy ✗ |
| **Nonlinear evolving** | R(t) = R₀ · f(E(t)) | Complex dynamics ✗ |

Where f can be sigmoid, tanh, or exponential decay — and E(t) is the flow energy at time t.

### Core Capabilities

- **Spectral gap tracking** over time for all sheaf models
- **Gap change rate**: dλ₁/dt — how fast is the gap evolving?
- **Phase transition detection**: when dλ₁/dt changes sign, something fundamental shifts
- **Stability maps**: for a given graph topology, map the (α, time) parameter space — where does the gap grow vs shrink?

## Build & Test

```bash
# Build
gcc -O2 -I include -o test_runner tests/test_evolving_sheaf.c src/evolving_sheaf.c -lm

# Run
./test_runner
```

## Quick Example

```c
#include "evolving_sheaf.h"

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

    es_gap_trajectory_free(traj);
    es_graph_free(g);
}
```

## Research Questions

This library is a tool for investigating open questions:

### 1. Is there a universal bound on gap decrease rate?

For linear evolving sheaves R(t) = R₀ + α·E(t), does |dλ₁/dt| have a universal upper bound in terms of α and the graph's Cheeger constant? Preliminary evidence suggests the bound is O(α·h(G)⁻¹).

### 2. Which topologies resist gap decrease?

**Hypothesis**: Expander graphs are robust — their high Cheeger constant means the spectral gap barely moves even under strong flow perturbation. Path graphs (h=1/n) are maximally fragile.

Early experiments support this: on 10-vertex graphs with α=1.0 and sinusoidal flow, the expander's relative gap change is significantly smaller than the cycle's.

### 3. Is there a "gap conservation law"?

When the gap decreases in one region of the spectrum, does it increase in another? If λ₁ drops, does λ₂ compensate? This would be a kind of spectral equipartition — a conservation law hidden in the Hodge structure.

### 4. Phase transitions as critical phenomena

The sign change in dλ₁/dt marks a phase transition in the sheaf's connectivity. What is the order parameter? Is it continuous or first-order? The sigmoid nonlinear model shows particularly rich transition behavior.

## Architecture

```
include/evolving_sheaf.h    — Public API (full header docs)
src/evolving_sheaf.c        — Implementation (Hodge Laplacian + eigenvalue solver)
tests/test_evolving_sheaf.c — 42 tests covering all models and edge cases
```

### Eigenvalue Solver

Uses Householder tridiagonalization → QL-iteration with implicit shifts. No external dependencies. All eigenvalues computed to ~1e-10 precision.

## Test Coverage

| Category | Tests | What's Verified |
|----------|-------|-----------------|
| Graph construction | 5 | Cycle, path, complete, expander |
| Static sheaf | 7 | Gap constant, known values, R₀ scaling |
| Linear evolving | 6 | Gap variation, α effects, pulse flow |
| Nonlinear evolving | 7 | Sigmoid, tanh, exp-decay |
| Phase transitions | 4 | Detection, static=none, rates |
| Stability maps | 2 | Grid computation, α scaling |
| Topology comparison | 3 | Expander robustness, path fragility |
| Restriction evaluation | 3 | Model-specific formulas |
| Edge cases | 5 | Single edge, large graphs, high α |
| Research observations | 1 | Gap sum experiment |
| **Total** | **43** | |

## The Philosophical Point

Most work on sheaf cohomology assumes static structure. The theorem is clean: spectral gap is an invariant. But nature isn't static — flows move, energy propagates, connectivity breathes.

When we let the sheaf evolve, the spectral gap becomes a *signal* rather than a constant. Its decrease tells us the system is losing coherence. Phase transitions mark critical moments. And different topologies respond differently — some resist, some don't.

**The failure of the theorem IS the research program.**

## License

MIT
