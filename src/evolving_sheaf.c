/*
 * evolving-sheaf-c: core implementation
 *
 * Builds the sheaf Hodge Laplacian L = D†D where D is the coboundary
 * operator whose entries depend on restriction maps. Restriction maps
 * may be static (theorem holds) or dynamic (gap decreases → research!).
 *
 * Eigenvalue computation uses Jacobi rotation method for symmetric matrices.
 */

#include "evolving_sheaf.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ═══════════════════════════════════════════
 *  Robust eigenvalue solver (Jacobi rotation)
 * ═══════════════════════════════════════════ */

static void jacobi_eigenvalues(double *A, int n, double *evals)
{
    /* Copy A since we destroy it */
    double *M = (double *)malloc(n * n * sizeof(double));
    memcpy(M, A, n * n * sizeof(double));

    for (int sweep = 0; sweep < 100 * n; sweep++) {
        /* Find largest off-diagonal element */
        double max_off = 0.0;
        int p = 0, q = 1;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (fabs(M[i * n + j]) > max_off) {
                    max_off = fabs(M[i * n + j]);
                    p = i; q = j;
                }

        if (max_off < 1e-14 * n) break;

        /* Compute rotation angle */
        double app = M[p * n + p], aqq = M[q * n + q], apq = M[p * n + q];
        double theta;
        if (fabs(app - aqq) < 1e-30)
            theta = M_PI / 4.0;
        else
            theta = 0.5 * atan2(2.0 * apq, app - aqq);

        double c = cos(theta), s = sin(theta);

        /* Apply rotation: M' = G^T M G */
        for (int i = 0; i < n; i++) {
            if (i == p || i == q) continue;
            double mip = M[i * n + p], miq = M[i * n + q];
            M[i * n + p] = c * mip + s * miq;
            M[p * n + i] = M[i * n + p];
            M[i * n + q] = -s * mip + c * miq;
            M[q * n + i] = M[i * n + q];
        }

        double new_pp = c*c*app + 2*s*c*apq + s*s*aqq;
        double new_qq = s*s*app - 2*s*c*apq + c*c*aqq;
        double new_pq = 0.0; /* should be zero by construction */

        M[p * n + p] = new_pp;
        M[q * n + q] = new_qq;
        M[p * n + q] = new_pq;
        M[q * n + p] = new_pq;
    }

    /* Extract diagonal */
    for (int i = 0; i < n; i++) evals[i] = M[i * n + i];

    /* Sort ascending */
    for (int i = 0; i < n - 1; i++)
        for (int j = i + 1; j < n; j++)
            if (evals[j] < evals[i]) { double t = evals[i]; evals[i] = evals[j]; evals[j] = t; }

    free(M);
}

/* ═══════════════════════════════════════════
 *  Graph
 * ═══════════════════════════════════════════ */

es_graph *es_graph_create(int nv, int ne, const es_edge *edges)
{
    es_graph *g = (es_graph *)calloc(1, sizeof(es_graph));
    g->n_vertices = nv;
    g->n_edges = ne;
    g->edges = (es_edge *)malloc(ne * sizeof(es_edge));
    memcpy(g->edges, edges, ne * sizeof(es_edge));
    return g;
}

void es_graph_free(es_graph *g)
{
    if (!g) return;
    free(g->edges);
    free(g->_adj_cap);
    free(g);
}

es_graph *es_graph_cycle(int n)
{
    es_edge *e = (es_edge *)malloc(n * sizeof(es_edge));
    for (int i = 0; i < n; i++) {
        e[i].src = i;
        e[i].dst = (i + 1) % n;
        e[i].weight = 1.0;
    }
    es_graph *g = es_graph_create(n, n, e);
    free(e);
    return g;
}

es_graph *es_graph_path(int n)
{
    if (n < 2) return es_graph_cycle(n);
    es_edge *e = (es_edge *)malloc((n - 1) * sizeof(es_edge));
    for (int i = 0; i < n - 1; i++) {
        e[i].src = i;
        e[i].dst = i + 1;
        e[i].weight = 1.0;
    }
    es_graph *g = es_graph_create(n, n - 1, e);
    free(e);
    return g;
}

es_graph *es_graph_complete(int n)
{
    int ne = n * (n - 1) / 2;
    es_edge *e = (es_edge *)malloc(ne * sizeof(es_edge));
    int k = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            e[k].src = i; e[k].dst = j; e[k].weight = 1.0; k++;
        }
    es_graph *g = es_graph_create(n, ne, e);
    free(e);
    return g;
}

es_graph *es_graph_expander(int n)
{
    int max_e = n * 3;
    es_edge *e = (es_edge *)malloc(max_e * sizeof(es_edge));
    int ne = 0;
    for (int i = 0; i < n; i++) {
        int targets[3] = { (i + 1) % n, (i + 2) % n, (i + n / 3) % n };
        for (int d = 0; d < 3; d++) {
            int a = i, b = targets[d];
            if (a > b) { int t = a; a = b; b = t; }
            int dup = 0;
            for (int j = 0; j < ne; j++)
                if (e[j].src == a && e[j].dst == b) { dup = 1; break; }
            if (!dup) {
                e[ne].src = a; e[ne].dst = b; e[ne].weight = 1.0;
                ne++;
            }
        }
    }
    es_graph *g = es_graph_create(n, ne, e);
    free(e);
    return g;
}

/* ═══════════════════════════════════════════
 *  Flow energy functions
 * ═══════════════════════════════════════════ */

double es_flow_constant(int edge, double t, void *ctx)
{
    (void)edge; (void)t; (void)ctx;
    return 1.0;
}

double es_flow_sinusoidal(int edge, double t, void *ctx)
{
    (void)ctx;
    return 1.0 + 0.5 * sin(0.5 * t + edge * 0.3);
}

double es_flow_pulse(int edge, double t, void *ctx)
{
    (void)ctx;
    double period = 4.0;
    double phase = fmod(t, period);
    if (phase < 0) phase += period;
    double pw = 0.3 + 0.05 * (edge % 5);
    return (phase < pw) ? 3.0 : 0.5;
}

static double deterministic_rand(int seed)
{
    unsigned int s = (unsigned int)(seed * 2654435761u);
    s = ((s >> 16) ^ s) * 0x45d9f3bu;
    s = ((s >> 16) ^ s) * 0x45d9f3bu;
    s = (s >> 16) ^ s;
    return (double)(s & 0x7FFFFFFF) / (double)0x7FFFFFFF;
}

double es_flow_random_walk(int edge, double t, void *ctx)
{
    (void)ctx;
    return 0.5 + deterministic_rand((int)(t * 100) + edge * 137);
}

/* ═══════════════════════════════════════════
 *  Restriction map evaluation
 * ═══════════════════════════════════════════ */

double es_eval_nonlin(es_nonlin_fn fn, double x, double k)
{
    switch (fn) {
    case ES_SIGMOID:
        return 1.0 / (1.0 + exp(-k * x));
    case ES_TANH:
        return tanh(k * x);
    case ES_EXPDECAY:
        return exp(-k * x * x);
    default:
        return x;
    }
}

double es_eval_restriction(const es_sheaf_config *cfg, double flow_energy)
{
    switch (cfg->model) {
    case ES_STATIC:
        return cfg->R0;
    case ES_LINEAR:
        return cfg->R0 + cfg->alpha * flow_energy;
    case ES_NONLINEAR:
        return cfg->R0 * es_eval_nonlin(cfg->nonlin, flow_energy, cfg->nonlin_k);
    default:
        return cfg->R0;
    }
}

/* ═══════════════════════════════════════════
 *  Spectrum computation
 * ═══════════════════════════════════════════ */

es_spectrum *es_compute_spectrum(const es_graph *g,
                                 const es_sheaf_config *cfg,
                                 es_flow_fn flow,
                                 void *flow_ctx,
                                 double t)
{
    int n = g->n_vertices;
    double *L = (double *)calloc(n * n, sizeof(double));

    for (int e = 0; e < g->n_edges; e++) {
        int i = g->edges[e].src;
        int j = g->edges[e].dst;
        double E = flow(e, t, flow_ctx);
        double r = es_eval_restriction(cfg, E) * g->edges[e].weight;

        L[i * n + i] += r * r;
        L[j * n + j] += r * r;
        L[i * n + j] -= r * r;
        L[j * n + i] -= r * r;
    }

    double *evals = (double *)malloc(n * sizeof(double));
    jacobi_eigenvalues(L, n, evals);
    free(L);

    es_spectrum *s = (es_spectrum *)calloc(1, sizeof(es_spectrum));
    s->t = t;
    s->eigenvalues = evals;
    s->n_eigenvalues = n;
    s->max_eigenvalue = evals[n - 1];

    s->trace = 0.0;
    for (int i = 0; i < n; i++) s->trace += evals[i];

    /* Find λ₁: smallest *nonzero* eigenvalue */
    s->lambda1 = -1.0;
    s->lambda2 = -1.0;
    int found = 0;
    for (int i = 0; i < n; i++) {
        if (evals[i] > 1e-10) {
            if (!found) { s->lambda1 = evals[i]; found = 1; }
            else { s->lambda2 = evals[i]; break; }
        }
    }
    if (s->lambda2 < 0) s->lambda2 = s->lambda1;

    return s;
}

void es_spectrum_free(es_spectrum *s)
{
    if (!s) return;
    free(s->eigenvalues);
    free(s);
}

/* ═══════════════════════════════════════════
 *  Gap trajectory tracking
 * ═══════════════════════════════════════════ */

es_gap_trajectory *es_track_gap(const es_graph *g,
                                const es_sheaf_config *cfg,
                                es_flow_fn flow,
                                void *flow_ctx,
                                double t0, double t1, int n_steps)
{
    es_gap_trajectory *traj = (es_gap_trajectory *)calloc(1, sizeof(es_gap_trajectory));
    traj->n_points = n_steps + 1;
    traj->points = (es_gap_point *)calloc(traj->n_points, sizeof(es_gap_point));
    traj->min_gap = DBL_MAX;
    traj->max_gap = 0.0;

    double dt = (n_steps > 0) ? (t1 - t0) / n_steps : 0.0;
    double prev_gap = -1.0;
    double prev_rate = 0.0;
    int n_trans = 0;

    for (int i = 0; i <= n_steps; i++) {
        double t = t0 + i * dt;
        es_spectrum *sp = es_compute_spectrum(g, cfg, flow, flow_ctx, t);

        es_gap_point *pt = &traj->points[i];
        pt->t = t;
        pt->gap = sp->lambda1;

        if (pt->gap < traj->min_gap) traj->min_gap = pt->gap;
        if (pt->gap > traj->max_gap) traj->max_gap = pt->gap;

        if (i > 0) {
            pt->gap_rate = (pt->gap - prev_gap) / dt;
            if (i > 1 && prev_rate * pt->gap_rate < -1e-12) {
                pt->phase_transition = 1;
                n_trans++;
            }
        }

        prev_gap = pt->gap;
        if (i > 0) prev_rate = pt->gap_rate;

        es_spectrum_free(sp);
    }

    traj->n_transitions = n_trans;
    traj->total_change = fabs(traj->max_gap - traj->min_gap);

    return traj;
}

void es_gap_trajectory_free(es_gap_trajectory *t)
{
    if (!t) return;
    free(t->points);
    free(t);
}

/* ═══════════════════════════════════════════
 *  Stability map
 * ═══════════════════════════════════════════ */

int es_stability_map(const es_graph *g,
                     const es_sheaf_config *base_cfg,
                     es_flow_fn flow, void *flow_ctx,
                     const double *alphas, int n_alphas,
                     const double *times, int n_times,
                     double *gap_values)
{
    for (int i = 0; i < n_alphas; i++) {
        es_sheaf_config cfg = *base_cfg;
        cfg.alpha = alphas[i];
        for (int j = 0; j < n_times; j++) {
            es_spectrum *sp = es_compute_spectrum(g, &cfg, flow, flow_ctx, times[j]);
            gap_values[i * n_times + j] = sp->lambda1;
            es_spectrum_free(sp);
        }
    }
    return 0;
}
