#ifndef EVOLVING_SHEAF_H
#define EVOLVING_SHEAF_H

#include <stddef.h>

/*
 * evolving-sheaf-c: Spectral Gap Dynamics in Evolving Cellular Sheaves
 *
 * Studies how the Hodge Laplacian's spectral gap (λ₁) changes when
 * restriction maps become flow-dependent. Static sheaves preserve the
 * gap; dynamic ones reveal phase transitions.
 */

/* ── Graph topology ── */

typedef struct {
    int src, dst;       /* directed edge */
    double weight;      /* edge weight (default 1.0) */
} es_edge;

typedef struct {
    int n_vertices;
    int n_edges;
    es_edge *edges;
    /* adjacency built internally */
    int *_adj_cap;
} es_graph;

/* Build a graph. edges[] is copied. */
es_graph *es_graph_create(int n_vertices, int n_edges, const es_edge *edges);
void      es_graph_free(es_graph *g);

/* Preset topologies */
es_graph *es_graph_cycle(int n);      /* cycle graph C_n */
es_graph *es_graph_complete(int n);   /* K_n */
es_graph *es_graph_path(int n);       /* P_n */
es_graph *es_graph_expander(int n);   /* 3-regular expander-like (Paley-ish) */

/* ── Sheaf model ── */

typedef enum {
    ES_STATIC    = 0,  /* R₀ constant — theorem holds */
    ES_LINEAR    = 1,  /* R(t) = R₀ + α·E(t) */
    ES_NONLINEAR = 2   /* R(t) = R₀ · f(E(t)), f = sigmoid/tanh/etc. */
} es_model;

typedef enum {
    ES_SIGMOID = 0,
    ES_TANH    = 1,
    ES_EXPDECAY = 2
} es_nonlin_fn;

typedef struct {
    es_model     model;
    double       R0;           /* base restriction map magnitude */
    double       alpha;        /* linear scaling coefficient */
    es_nonlin_fn nonlin;       /* nonlinear function selector */
    double       nonlin_k;     /* steepness parameter for nonlinear */
} es_sheaf_config;

/* ── Flow energy ── */

typedef double (*es_flow_fn)(int edge_idx, double t, void *ctx);

/* Built-in flow generators */
double es_flow_constant(int edge, double t, void *ctx);
double es_flow_sinusoidal(int edge, double t, void *ctx);
double es_flow_pulse(int edge, double t, void *ctx);
double es_flow_random_walk(int edge, double t, void *ctx);

/* ── Spectral analysis ── */

typedef struct {
    double t;              /* time point */
    double lambda1;        /* smallest nonzero eigenvalue (spectral gap) */
    double lambda2;        /* second nonzero eigenvalue */
    double max_eigenvalue; /* largest eigenvalue */
    double trace;          /* trace of Laplacian */
    double *eigenvalues;   /* all eigenvalues (n_vertices) */
    int    n_eigenvalues;
} es_spectrum;

void es_spectrum_free(es_spectrum *s);

typedef struct {
    double t;
    double gap;
    double gap_rate;       /* dλ₁/dt (numerical) */
    int    phase_transition; /* 1 if sign of dλ₁/dt changed */
} es_gap_point;

typedef struct {
    int           n_points;
    es_gap_point *points;
    double        min_gap;       /* minimum spectral gap observed */
    double        max_gap;       /* maximum spectral gap observed */
    double        total_change;  /* |Δλ₁| over trajectory */
    int           n_transitions; /* phase transition count */
} es_gap_trajectory;

void es_gap_trajectory_free(es_gap_trajectory *t);

/* ── Core computation ── */

/*
 * Compute the Hodge Laplacian for the sheaf at time t.
 * Returns spectrum (caller frees via es_spectrum_free).
 */
es_spectrum *es_compute_spectrum(const es_graph *g,
                                 const es_sheaf_config *cfg,
                                 es_flow_fn flow,
                                 void *flow_ctx,
                                 double t);

/*
 * Track spectral gap over a time interval [t0, t1] with n_steps steps.
 * Returns trajectory (caller frees).
 */
es_gap_trajectory *es_track_gap(const es_graph *g,
                                const es_sheaf_config *cfg,
                                es_flow_fn flow,
                                void *flow_ctx,
                                double t0, double t1, int n_steps);

/*
 * Stability map: compute gap over (alpha, time) grid.
 * Outputs: gap_values[i * n_times + j] = gap at (alphas[i], times[j])
 * Returns 0 on success.
 */
int es_stability_map(const es_graph *g,
                     const es_sheaf_config *base_cfg,
                     es_flow_fn flow, void *flow_ctx,
                     const double *alphas, int n_alphas,
                     const double *times, int n_times,
                     double *gap_values);

/* ── Restriction map evaluation (exposed for testing) ── */

double es_eval_restriction(const es_sheaf_config *cfg, double flow_energy);
double es_eval_nonlin(es_nonlin_fn fn, double x, double k);

/* ── Version ── */

#define EVOLVING_SHEAF_VERSION "0.1.0"
#define EVOLVING_SHEAF_VERSION_MAJOR 0
#define EVOLVING_SHEAF_VERSION_MINOR 1
#define EVOLVING_SHEAF_VERSION_PATCH 0

#endif /* EVOLVING_SHEAF_H */
