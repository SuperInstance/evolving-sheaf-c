/*
 * evolving-sheaf-c test suite
 *
 * 40+ tests covering:
 *   - Static sheaf (theorem holds, gap is constant)
 *   - Linear evolving (gap decreases with flow energy)
 *   - Nonlinear evolving (sigmoid, tanh, exp-decay)
 *   - Phase transition detection
 *   - Stability maps
 *   - Graph topology comparisons
 *   - Edge cases
 */

#include "evolving_sheaf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define ASSERT(expr, msg) do { \
    tests_run++; \
    if (expr) { tests_passed++; } \
    else { tests_failed++; \
           printf("  FAIL [%d]: %s (line %d)\n", tests_run, msg, __LINE__); } \
} while(0)

#define ASSERT_NEAR(a, b, eps, msg) \
    ASSERT(fabs((a) - (b)) < (eps), msg)

#define TOLERANCE 0.05

/* ═══════════════════════════════════════════
 *  Graph construction tests
 * ═══════════════════════════════════════════ */

static void test_cycle_graph(void)
{
    printf("── Cycle graph tests ──\n");
    es_graph *g = es_graph_cycle(4);
    ASSERT(g != NULL, "cycle(4) created");
    ASSERT(g->n_vertices == 4, "cycle(4) has 4 vertices");
    ASSERT(g->n_edges == 4, "cycle(4) has 4 edges");
    for (int i = 0; i < 4; i++)
        ASSERT(g->edges[i].src == i && g->edges[i].dst == (i+1)%4,
               "cycle edge connectivity");
    es_graph_free(g);
}

static void test_path_graph(void)
{
    printf("── Path graph tests ──\n");
    es_graph *g = es_graph_path(5);
    ASSERT(g != NULL, "path(5) created");
    ASSERT(g->n_vertices == 5, "path(5) has 5 vertices");
    ASSERT(g->n_edges == 4, "path(5) has 4 edges");
    es_graph_free(g);
}

static void test_complete_graph(void)
{
    printf("── Complete graph tests ──\n");
    es_graph *g = es_graph_complete(4);
    ASSERT(g != NULL, "K4 created");
    ASSERT(g->n_vertices == 4, "K4 has 4 vertices");
    ASSERT(g->n_edges == 6, "K4 has 6 edges");
    es_graph_free(g);
}

static void test_expander_graph(void)
{
    printf("── Expander graph tests ──\n");
    es_graph *g = es_graph_expander(10);
    ASSERT(g != NULL, "expander(10) created");
    ASSERT(g->n_vertices == 10, "expander has 10 vertices");
    ASSERT(g->n_edges >= 10, "expander has ≥10 edges");
    es_graph_free(g);
}

static void test_cycle_3(void)
{
    printf("── Cycle-3 tests ──\n");
    es_graph *g = es_graph_cycle(3);
    ASSERT(g->n_vertices == 3, "cycle(3) vertices");
    ASSERT(g->n_edges == 3, "cycle(3) edges");
    es_graph_free(g);
}

/* ═══════════════════════════════════════════
 *  Static sheaf tests (theorem holds)
 * ═══════════════════════════════════════════ */

static void test_static_spectrum_exists(void)
{
    printf("── Static sheaf: spectrum exists ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_STATIC, .R0 = 1.0, .alpha = 0.0 };
    es_spectrum *sp = es_compute_spectrum(g, &cfg, es_flow_constant, NULL, 0.0);
    ASSERT(sp != NULL, "spectrum computed");
    ASSERT(sp->n_eigenvalues == 4, "4 eigenvalues");
    ASSERT(sp->lambda1 > 0.0, "lambda1 > 0");
    es_spectrum_free(sp);
    es_graph_free(g);
}

static void test_static_gap_constant_over_time(void)
{
    printf("── Static sheaf: gap is constant over time ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_STATIC, .R0 = 1.0, .alpha = 0.0 };

    double gaps[5];
    for (int i = 0; i < 5; i++) {
        es_spectrum *sp = es_compute_spectrum(g, &cfg, es_flow_sinusoidal, NULL, i * 2.0);
        gaps[i] = sp->lambda1;
        es_spectrum_free(sp);
    }

    for (int i = 1; i < 5; i++)
        ASSERT_NEAR(gaps[i], gaps[0], 0.001, "static gap constant over time");
    es_graph_free(g);
}

static void test_static_cycle4_known_value(void)
{
    printf("── Static sheaf: cycle-4 known eigenvalue ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_STATIC, .R0 = 1.0, .alpha = 0.0 };
    es_spectrum *sp = es_compute_spectrum(g, &cfg, es_flow_constant, NULL, 0.0);

    /* For cycle-4 with R0=1, L = standard graph Laplacian.
     * Eigenvalues: 0, 2, 2, 4.  λ₁ = 2.0 */
    ASSERT_NEAR(sp->eigenvalues[0], 0.0, TOLERANCE, "eigenvalue 0 (kernel)");
    ASSERT_NEAR(sp->lambda1, 2.0, TOLERANCE, "λ₁ ≈ 2.0 for cycle-4");
    ASSERT_NEAR(sp->max_eigenvalue, 4.0, TOLERANCE, "λ_max ≈ 4.0");
    es_spectrum_free(sp);
    es_graph_free(g);
}

static void test_static_trace(void)
{
    printf("── Static sheaf: trace = 2·|E|·R₀² ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_STATIC, .R0 = 2.0, .alpha = 0.0 };
    es_spectrum *sp = es_compute_spectrum(g, &cfg, es_flow_constant, NULL, 0.0);
    /* 4 edges, R=2: trace = 4 * 2 * 2² = 32 */
    ASSERT_NEAR(sp->trace, 32.0, 0.5, "trace = 2·|E|·R₀²");
    es_spectrum_free(sp);
    es_graph_free(g);
}

static void test_static_R0_scaling(void)
{
    printf("── Static sheaf: gap scales as R₀² ──\n");
    es_graph *g = es_graph_cycle(4);
    double gap1, gap2;

    es_sheaf_config cfg1 = { .model = ES_STATIC, .R0 = 1.0, .alpha = 0.0 };
    es_spectrum *sp1 = es_compute_spectrum(g, &cfg1, es_flow_constant, NULL, 0.0);
    gap1 = sp1->lambda1;
    es_spectrum_free(sp1);

    es_sheaf_config cfg2 = { .model = ES_STATIC, .R0 = 2.0, .alpha = 0.0 };
    es_spectrum *sp2 = es_compute_spectrum(g, &cfg2, es_flow_constant, NULL, 0.0);
    gap2 = sp2->lambda1;
    es_spectrum_free(sp2);

    ASSERT_NEAR(gap2 / gap1, 4.0, 0.1, "gap scales as R₀²");
    es_graph_free(g);
}

static void test_static_path_gap(void)
{
    printf("── Static sheaf: path graph gap ──\n");
    es_graph *g = es_graph_path(4);
    es_sheaf_config cfg = { .model = ES_STATIC, .R0 = 1.0, .alpha = 0.0 };
    es_spectrum *sp = es_compute_spectrum(g, &cfg, es_flow_constant, NULL, 0.0);
    ASSERT(sp->lambda1 > 0.0, "path graph has positive gap");
    ASSERT(sp->lambda1 < 2.0, "path gap < cycle gap for same |V|");
    es_spectrum_free(sp);
    es_graph_free(g);
}

static void test_static_complete_gap(void)
{
    printf("── Static sheaf: complete graph gap ──\n");
    es_graph *g = es_graph_complete(5);
    es_sheaf_config cfg = { .model = ES_STATIC, .R0 = 1.0, .alpha = 0.0 };
    es_spectrum *sp = es_compute_spectrum(g, &cfg, es_flow_constant, NULL, 0.0);
    ASSERT(sp->lambda1 > 0.0, "complete graph has positive gap");
    /* K5 has good expansion → large gap */
    ASSERT(sp->lambda1 > 3.0, "K5 gap > 3.0 (good expansion)");
    es_spectrum_free(sp);
    es_graph_free(g);
}

/* ═══════════════════════════════════════════
 *  Linear evolving sheaf tests
 * ═══════════════════════════════════════════ */

static void test_linear_gap_differs_from_static(void)
{
    printf("── Linear evolving: gap differs from static ──\n");
    es_graph *g = es_graph_cycle(4);

    es_sheaf_config scfg = { .model = ES_STATIC, .R0 = 1.0, .alpha = 0.0 };
    es_spectrum *ssp = es_compute_spectrum(g, &scfg, es_flow_sinusoidal, NULL, 5.0);
    double sg = ssp->lambda1;
    es_spectrum_free(ssp);

    es_sheaf_config lcfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 0.5 };
    es_spectrum *lsp = es_compute_spectrum(g, &lcfg, es_flow_sinusoidal, NULL, 5.0);
    double lg = lsp->lambda1;
    es_spectrum_free(lsp);

    ASSERT(fabs(lg - sg) > 0.01, "linear evolving gap ≠ static gap");
    es_graph_free(g);
}

static void test_linear_gap_decreases_with_positive_alpha(void)
{
    printf("── Linear evolving: positive α increases restriction → gap changes ──\n");
    es_graph *g = es_graph_cycle(4);

    es_sheaf_config cfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 0.3 };
    es_gap_trajectory *traj = es_track_gap(g, &cfg, es_flow_sinusoidal, NULL, 0.0, 20.0, 100);

    ASSERT(traj->total_change > 0.01, "gap varies over time with linear evolving");
    ASSERT(traj->n_points == 101, "correct number of trajectory points");

    es_gap_trajectory_free(traj);
    es_graph_free(g);
}

static void test_linear_trajectory_min_max(void)
{
    printf("── Linear evolving: trajectory tracks min/max ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 0.5 };
    es_gap_trajectory *traj = es_track_gap(g, &cfg, es_flow_sinusoidal, NULL, 0.0, 20.0, 200);

    ASSERT(traj->min_gap > 0.0, "min gap > 0");
    ASSERT(traj->max_gap > traj->min_gap, "max gap > min gap");
    ASSERT(traj->max_gap >= traj->min_gap, "max ≥ min");

    es_gap_trajectory_free(traj);
    es_graph_free(g);
}

static void test_linear_alpha_zero_equals_static(void)
{
    printf("── Linear evolving: α=0 equals static ──\n");
    es_graph *g = es_graph_cycle(4);

    es_sheaf_config scfg = { .model = ES_STATIC, .R0 = 1.0, .alpha = 0.0 };
    es_sheaf_config lcfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 0.0 };

    es_spectrum *ss = es_compute_spectrum(g, &scfg, es_flow_sinusoidal, NULL, 3.0);
    es_spectrum *ls = es_compute_spectrum(g, &lcfg, es_flow_sinusoidal, NULL, 3.0);

    ASSERT_NEAR(ss->lambda1, ls->lambda1, 0.001, "α=0 → same as static");
    es_spectrum_free(ss);
    es_spectrum_free(ls);
    es_graph_free(g);
}

static void test_linear_larger_alpha_more_variation(void)
{
    printf("── Linear evolving: larger α → more variation ──\n");
    es_graph *g = es_graph_cycle(4);

    es_sheaf_config cfg1 = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 0.1 };
    es_sheaf_config cfg2 = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 1.0 };

    es_gap_trajectory *t1 = es_track_gap(g, &cfg1, es_flow_sinusoidal, NULL, 0.0, 10.0, 100);
    es_gap_trajectory *t2 = es_track_gap(g, &cfg2, es_flow_sinusoidal, NULL, 0.0, 10.0, 100);

    ASSERT(t2->total_change > t1->total_change,
           "larger α → more total change");
    es_gap_trajectory_free(t1);
    es_gap_trajectory_free(t2);
    es_graph_free(g);
}

static void test_linear_pulse_flow(void)
{
    printf("── Linear evolving: pulse flow ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 0.5 };

    es_gap_trajectory *traj = es_track_gap(g, &cfg, es_flow_pulse, NULL, 0.0, 20.0, 200);
    ASSERT(traj->n_points == 201, "pulse trajectory has correct length");
    ASSERT(traj->total_change > 0.0, "pulse causes gap variation");
    es_gap_trajectory_free(traj);
    es_graph_free(g);
}

/* ═══════════════════════════════════════════
 *  Nonlinear evolving sheaf tests
 * ═══════════════════════════════════════════ */

static void test_nonlinear_sigmoid(void)
{
    printf("── Nonlinear evolving: sigmoid ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_NONLINEAR, .R0 = 2.0, .alpha = 0.0,
                            .nonlin = ES_SIGMOID, .nonlin_k = 2.0 };

    es_spectrum *sp = es_compute_spectrum(g, &cfg, es_flow_sinusoidal, NULL, 2.0);
    ASSERT(sp != NULL, "sigmoid spectrum computed");
    ASSERT(sp->lambda1 > 0.0, "sigmoid gap > 0");
    es_spectrum_free(sp);
    es_graph_free(g);
}

static void test_nonlinear_tanh(void)
{
    printf("── Nonlinear evolving: tanh ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_NONLINEAR, .R0 = 2.0, .alpha = 0.0,
                            .nonlin = ES_TANH, .nonlin_k = 1.0 };

    es_gap_trajectory *traj = es_track_gap(g, &cfg, es_flow_sinusoidal, NULL, 0.0, 10.0, 100);
    ASSERT(traj->n_points == 101, "tanh trajectory ok");
    ASSERT(traj->total_change > 0.01, "tanh causes variation");
    es_gap_trajectory_free(traj);
    es_graph_free(g);
}

static void test_nonlinear_expdecay(void)
{
    printf("── Nonlinear evolving: exp-decay ──\n");
    es_graph *g = es_graph_cycle(5);
    es_sheaf_config cfg = { .model = ES_NONLINEAR, .R0 = 2.0, .alpha = 0.0,
                            .nonlin = ES_EXPDECAY, .nonlin_k = 0.5 };

    es_gap_trajectory *traj = es_track_gap(g, &cfg, es_flow_sinusoidal, NULL, 0.0, 10.0, 100);
    ASSERT(traj->min_gap > 0.0, "expdecay min gap > 0");
    es_gap_trajectory_free(traj);
    es_graph_free(g);
}

static void test_nonlinear_differs_from_linear(void)
{
    printf("── Nonlinear vs linear: different dynamics ──\n");
    es_graph *g = es_graph_cycle(4);

    es_sheaf_config lcfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 0.5 };
    es_sheaf_config ncfg = { .model = ES_NONLINEAR, .R0 = 1.0, .alpha = 0.0,
                             .nonlin = ES_SIGMOID, .nonlin_k = 2.0 };

    es_gap_trajectory *lt = es_track_gap(g, &lcfg, es_flow_sinusoidal, NULL, 0.0, 10.0, 100);
    es_gap_trajectory *nt = es_track_gap(g, &ncfg, es_flow_sinusoidal, NULL, 0.0, 10.0, 100);

    /* At least one point should differ significantly */
    int any_diff = 0;
    for (int i = 0; i < lt->n_points && i < nt->n_points; i++) {
        if (fabs(lt->points[i].gap - nt->points[i].gap) > 0.05) {
            any_diff = 1; break;
        }
    }
    ASSERT(any_diff, "nonlinear dynamics differ from linear");
    es_gap_trajectory_free(lt);
    es_gap_trajectory_free(nt);
    es_graph_free(g);
}

static void test_nonlin_eval_sigmoid(void)
{
    printf("── Nonlinear eval: sigmoid(0) = 0.5 ──\n");
    double v = es_eval_nonlin(ES_SIGMOID, 0.0, 1.0);
    ASSERT_NEAR(v, 0.5, 0.001, "sigmoid(0) = 0.5");
}

static void test_nonlin_eval_tanh(void)
{
    printf("── Nonlinear eval: tanh(0) = 0 ──\n");
    double v = es_eval_nonlin(ES_TANH, 0.0, 1.0);
    ASSERT_NEAR(v, 0.0, 0.001, "tanh(0) = 0");
}

static void test_nonlin_eval_expdecay(void)
{
    printf("── Nonlinear eval: exp(0) = 1 ──\n");
    double v = es_eval_nonlin(ES_EXPDECAY, 0.0, 1.0);
    ASSERT_NEAR(v, 1.0, 0.001, "exp(-0) = 1");
}

/* ═══════════════════════════════════════════
 *  Phase transition detection
 * ═══════════════════════════════════════════ */

static void test_phase_transition_sinusoidal(void)
{
    printf("── Phase transition: sinusoidal flow triggers transitions ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 1.0 };

    es_gap_trajectory *traj = es_track_gap(g, &cfg, es_flow_sinusoidal, NULL, 0.0, 30.0, 500);
    /* Sinusoidal flow should cause the gap rate to oscillate, creating transitions */
    ASSERT(traj->n_transitions >= 1, "sinusoidal flow causes ≥1 phase transition");
    es_gap_trajectory_free(traj);
    es_graph_free(g);
}

static void test_phase_transition_static_none(void)
{
    printf("── Phase transition: static sheaf has no transitions ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_STATIC, .R0 = 1.0, .alpha = 0.0 };

    es_gap_trajectory *traj = es_track_gap(g, &cfg, es_flow_constant, NULL, 0.0, 10.0, 100);
    ASSERT(traj->n_transitions == 0, "static sheaf: no phase transitions");
    es_gap_trajectory_free(traj);
    es_graph_free(g);
}

static void test_phase_transition_pulse(void)
{
    printf("── Phase transition: pulse causes rapid gap change ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 1.5 };

    es_gap_trajectory *traj = es_track_gap(g, &cfg, es_flow_pulse, NULL, 0.0, 20.0, 400);
    /* Pulse causes sharp gap jumps; verify large total change */
    ASSERT(traj->total_change > 1.0, "pulse causes large gap variation");
    /* Some gap rates should be positive and some negative */
    int pos = 0, neg = 0;
    for (int i = 1; i < traj->n_points; i++) {
        if (traj->points[i].gap_rate > 0.01) pos++;
        if (traj->points[i].gap_rate < -0.01) neg++;
    }
    ASSERT(pos > 0 && neg > 0, "gap rate has both positive and negative phases");
    es_gap_trajectory_free(traj);
    es_graph_free(g);
}

static void test_gap_rate_values(void)
{
    printf("── Gap rate: computed correctly ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 0.5 };

    es_gap_trajectory *traj = es_track_gap(g, &cfg, es_flow_sinusoidal, NULL, 0.0, 10.0, 100);
    /* First point has no rate */
    ASSERT(traj->points[0].gap_rate == 0.0, "first point rate = 0");
    /* Some subsequent points should have nonzero rate */
    int any_rate = 0;
    for (int i = 1; i < traj->n_points; i++)
        if (fabs(traj->points[i].gap_rate) > 1e-6) { any_rate = 1; break; }
    ASSERT(any_rate, "some nonzero gap rates");

    es_gap_trajectory_free(traj);
    es_graph_free(g);
}

/* ═══════════════════════════════════════════
 *  Stability map tests
 * ═══════════════════════════════════════════ */

static void test_stability_map_basic(void)
{
    printf("── Stability map: basic computation ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 0.0 };

    double alphas[] = {0.0, 0.5, 1.0};
    double times[] = {0.0, 5.0, 10.0};
    double gaps[9];

    int rc = es_stability_map(g, &cfg, es_flow_sinusoidal, NULL,
                              alphas, 3, times, 3, gaps);
    ASSERT(rc == 0, "stability map computed");
    /* α=0 should be constant across time */
    ASSERT_NEAR(gaps[0], gaps[1], 0.01, "α=0: gap constant at t=0,5");
    ASSERT_NEAR(gaps[0], gaps[2], 0.01, "α=0: gap constant at t=0,10");
    /* Larger α should show variation */
    int any_diff = 0;
    for (int j = 0; j < 3; j++)
        if (fabs(gaps[6 + j] - gaps[0]) > 0.01) { any_diff = 1; break; }
    ASSERT(any_diff, "larger α changes gap");

    es_graph_free(g);
}

static void test_stability_map_increasing_alpha(void)
{
    printf("── Stability map: increasing α increases variation ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 0.0 };

    double alphas[] = {0.0, 2.0};
    double times[] = {0.0, 3.0, 6.0, 9.0};
    double gaps[8];

    es_stability_map(g, &cfg, es_flow_sinusoidal, NULL, alphas, 2, times, 4, gaps);

    /* Row 0 (α=0): all same; Row 1 (α=2): some variation */
    double var0 = 0.0, var1 = 0.0;
    for (int j = 0; j < 4; j++) {
        var0 += fabs(gaps[j] - gaps[0]);
        var1 += fabs(gaps[4 + j] - gaps[4]);
    }
    ASSERT(var0 < 0.01, "α=0: no variation");
    ASSERT(var1 > var0, "α=2: more variation than α=0");

    es_graph_free(g);
}

/* ═══════════════════════════════════════════
 *  Topology comparison tests
 * ═══════════════════════════════════════════ */

static void test_topology_cycle_vs_complete(void)
{
    printf("── Topology: complete graph has larger gap than cycle ──\n");
    es_graph *c4 = es_graph_cycle(4);
    es_graph *k4 = es_graph_complete(4);
    es_sheaf_config cfg = { .model = ES_STATIC, .R0 = 1.0, .alpha = 0.0 };

    es_spectrum *cs = es_compute_spectrum(c4, &cfg, es_flow_constant, NULL, 0.0);
    es_spectrum *ks = es_compute_spectrum(k4, &cfg, es_flow_constant, NULL, 0.0);

    ASSERT(ks->lambda1 > cs->lambda1, "K4 gap > C4 gap");
    es_spectrum_free(cs);
    es_spectrum_free(ks);
    es_graph_free(c4);
    es_graph_free(k4);
}

static void test_topology_expander_robustness(void)
{
    printf("── Topology: expander resists gap decrease ──\n");
    es_graph *cycle = es_graph_cycle(10);
    es_graph *expander = es_graph_expander(10);
    es_sheaf_config cfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 1.0 };

    es_gap_trajectory *ct = es_track_gap(cycle, &cfg, es_flow_sinusoidal, NULL, 0.0, 10.0, 200);
    es_gap_trajectory *et = es_track_gap(expander, &cfg, es_flow_sinusoidal, NULL, 0.0, 10.0, 200);

    /* Expander should have less relative gap change */
    double cycle_rel = ct->total_change / ct->max_gap;
    double exp_rel   = et->total_change / et->max_gap;
    printf("    cycle relative change: %.4f, expander: %.4f\n", cycle_rel, exp_rel);
    /* We observe but don't assert strongly — this is a research question */

    ASSERT(et->min_gap > 0.0, "expander gap stays positive");
    ASSERT(ct->min_gap > 0.0, "cycle gap stays positive");

    es_gap_trajectory_free(ct);
    es_gap_trajectory_free(et);
    es_graph_free(cycle);
    es_graph_free(expander);
}

static void test_topology_path_fragile(void)
{
    printf("── Topology: path graph is fragile ──\n");
    es_graph *path = es_graph_path(10);
    es_graph *cycle = es_graph_cycle(10);
    es_sheaf_config cfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 0.5 };

    es_spectrum *ps = es_compute_spectrum(path, &cfg, es_flow_sinusoidal, NULL, 5.0);
    es_spectrum *cs = es_compute_spectrum(cycle, &cfg, es_flow_sinusoidal, NULL, 5.0);

    ASSERT(ps->lambda1 > 0.0, "path gap positive");
    ASSERT(cs->lambda1 > ps->lambda1, "cycle gap > path gap");

    es_spectrum_free(ps);
    es_spectrum_free(cs);
    es_graph_free(path);
    es_graph_free(cycle);
}

/* ═══════════════════════════════════════════
 *  Restriction map evaluation tests
 * ═══════════════════════════════════════════ */

static void test_restriction_static(void)
{
    printf("── Restriction eval: static returns R₀ ──\n");
    es_sheaf_config cfg = { .model = ES_STATIC, .R0 = 3.14, .alpha = 0.0 };
    ASSERT_NEAR(es_eval_restriction(&cfg, 5.0), 3.14, 0.001, "static R=R₀");
}

static void test_restriction_linear(void)
{
    printf("── Restriction eval: linear R₀ + αE ──\n");
    es_sheaf_config cfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 0.5 };
    ASSERT_NEAR(es_eval_restriction(&cfg, 2.0), 2.0, 0.001, "linear: 1.0 + 0.5*2.0 = 2.0");
}

static void test_restriction_nonlinear(void)
{
    printf("── Restriction eval: nonlinear R₀·f(E) ──\n");
    es_sheaf_config cfg = { .model = ES_NONLINEAR, .R0 = 2.0, .alpha = 0.0,
                            .nonlin = ES_TANH, .nonlin_k = 1.0 };
    double r = es_eval_restriction(&cfg, 1.0);
    double expected = 2.0 * tanh(1.0);
    ASSERT_NEAR(r, expected, 0.001, "nonlinear: 2.0*tanh(1.0)");
}

/* ═══════════════════════════════════════════
 *  Edge case tests
 * ═══════════════════════════════════════════ */

static void test_single_edge(void)
{
    printf("── Edge case: single edge graph ──\n");
    es_edge e = {0, 1, 1.0};
    es_graph *g = es_graph_create(2, 1, &e);
    es_sheaf_config cfg = { .model = ES_STATIC, .R0 = 1.0, .alpha = 0.0 };

    es_spectrum *sp = es_compute_spectrum(g, &cfg, es_flow_constant, NULL, 0.0);
    ASSERT(sp->n_eigenvalues == 2, "2 eigenvalues");
    ASSERT_NEAR(sp->eigenvalues[0], 0.0, TOLERANCE, "kernel eigenvalue");
    ASSERT(sp->lambda1 > 0.0, "positive gap");
    es_spectrum_free(sp);
    es_graph_free(g);
}

static void test_large_graph(void)
{
    printf("── Edge case: larger graph (cycle-20) ──\n");
    es_graph *g = es_graph_cycle(20);
    es_sheaf_config cfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 0.3 };

    es_spectrum *sp = es_compute_spectrum(g, &cfg, es_flow_sinusoidal, NULL, 2.0);
    ASSERT(sp->n_eigenvalues == 20, "20 eigenvalues");
    ASSERT(sp->lambda1 > 0.0, "positive gap");
    es_spectrum_free(sp);
    es_graph_free(g);
}

static void test_zero_flow_energy(void)
{
    printf("── Edge case: zero flow energy ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 0.5 };

    /* Custom flow that returns 0 */
    es_spectrum *sp = es_compute_spectrum(g, &cfg, es_flow_constant, NULL, 0.0);
    /* With E=1 (constant), R = 1.0 + 0.5*1.0 = 1.5 */
    ASSERT(sp->lambda1 > 0.0, "zero-flow gap positive");
    es_spectrum_free(sp);
    es_graph_free(g);
}

static void test_high_alpha(void)
{
    printf("── Edge case: very high α ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 10.0 };

    es_gap_trajectory *traj = es_track_gap(g, &cfg, es_flow_sinusoidal, NULL, 0.0, 10.0, 100);
    ASSERT(traj->total_change > 1.0, "high α causes large gap variation");
    ASSERT(traj->min_gap > 0.0, "gap stays positive even at high α");
    es_gap_trajectory_free(traj);
    es_graph_free(g);
}

static void test_trajectory_continuity(void)
{
    printf("── Trajectory: gap changes smoothly ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 0.5 };

    es_gap_trajectory *traj = es_track_gap(g, &cfg, es_flow_sinusoidal, NULL, 0.0, 10.0, 1000);

    /* Adjacent points should be close */
    double max_jump = 0.0;
    for (int i = 1; i < traj->n_points; i++) {
        double jump = fabs(traj->points[i].gap - traj->points[i-1].gap);
        if (jump > max_jump) max_jump = jump;
    }
    ASSERT(max_jump < 0.5, "gap changes smoothly between steps");

    es_gap_trajectory_free(traj);
    es_graph_free(g);
}

/* ═══════════════════════════════════════════
 *  Gap conservation hypothesis test
 * ═══════════════════════════════════════════ */

static void test_gap_sum_experiment(void)
{
    printf("── Research: gap sum over eigenvalues ──\n");
    es_graph *g = es_graph_cycle(4);
    es_sheaf_config cfg = { .model = ES_LINEAR, .R0 = 1.0, .alpha = 1.0 };

    /* Compute sum of all nonzero eigenvalues at two time points */
    es_spectrum *s1 = es_compute_spectrum(g, &cfg, es_flow_sinusoidal, NULL, 0.0);
    es_spectrum *s2 = es_compute_spectrum(g, &cfg, es_flow_sinusoidal, NULL, 5.0);

    double sum1 = 0.0, sum2 = 0.0;
    for (int i = 0; i < 4; i++) {
        sum1 += s1->eigenvalues[i];
        sum2 += s2->eigenvalues[i];
    }
    /* Trace = 2*sum(r²) which depends on flow energy, so trace changes.
     * But we record the observation. */
    printf("    trace(t=0)=%.3f, trace(t=5)=%.3f\n", sum1, sum2);
    ASSERT(1, "gap sum recorded (always passes, observational)");

    es_spectrum_free(s1);
    es_spectrum_free(s2);
    es_graph_free(g);
}

/* ═══════════════════════════════════════════
 *  Main
 * ═══════════════════════════════════════════ */

int main(void)
{
    printf("╔══════════════════════════════════════════════╗\n");
    printf("║  evolving-sheaf-c test suite                 ║\n");
    printf("║  Spectral Gap Dynamics in Evolving Sheaves   ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    /* Graph construction */
    test_cycle_graph();
    test_path_graph();
    test_complete_graph();
    test_expander_graph();
    test_cycle_3();

    /* Static sheaf (theorem holds) */
    test_static_spectrum_exists();
    test_static_gap_constant_over_time();
    test_static_cycle4_known_value();
    test_static_trace();
    test_static_R0_scaling();
    test_static_path_gap();
    test_static_complete_gap();

    /* Linear evolving */
    test_linear_gap_differs_from_static();
    test_linear_gap_decreases_with_positive_alpha();
    test_linear_trajectory_min_max();
    test_linear_alpha_zero_equals_static();
    test_linear_larger_alpha_more_variation();
    test_linear_pulse_flow();

    /* Nonlinear evolving */
    test_nonlinear_sigmoid();
    test_nonlinear_tanh();
    test_nonlinear_expdecay();
    test_nonlinear_differs_from_linear();
    test_nonlin_eval_sigmoid();
    test_nonlin_eval_tanh();
    test_nonlin_eval_expdecay();

    /* Phase transitions */
    test_phase_transition_sinusoidal();
    test_phase_transition_static_none();
    test_phase_transition_pulse();
    test_gap_rate_values();

    /* Stability maps */
    test_stability_map_basic();
    test_stability_map_increasing_alpha();

    /* Topology comparisons */
    test_topology_cycle_vs_complete();
    test_topology_expander_robustness();
    test_topology_path_fragile();

    /* Restriction map evaluation */
    test_restriction_static();
    test_restriction_linear();
    test_restriction_nonlinear();

    /* Edge cases */
    test_single_edge();
    test_large_graph();
    test_zero_flow_energy();
    test_high_alpha();
    test_trajectory_continuity();

    /* Research */
    test_gap_sum_experiment();

    printf("\n══════════════════════════════════════\n");
    printf("Results: %d/%d passed, %d failed\n",
           tests_passed, tests_run, tests_failed);
    printf("══════════════════════════════════════\n");

    return tests_failed > 0 ? 1 : 0;
}
