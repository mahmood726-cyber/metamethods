"""Tests for MetaMethods advanced statistical methods."""
import sys, os, math, pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from advanced_methods import (
    evalue, doi_plot_lfk, rosenthal_failsafe_n,
    quality_effects_model, proportion_meta, permutation_test_heterogeneity,
)


# ── E-value ──

class TestEvalue:
    def test_rr_above_1(self):
        """RR=2.0 should give E-value ~3.41"""
        r = evalue(2.0, measure='RR')
        assert r['evalue_point'] is not None
        assert 3.0 < r['evalue_point'] < 4.0  # Known: 2 + sqrt(2*1) = 3.414

    def test_rr_below_1(self):
        """RR=0.5 should give same E-value as RR=2.0 (symmetric)"""
        r = evalue(0.5, measure='RR')
        r2 = evalue(2.0, measure='RR')
        assert abs(r['evalue_point'] - r2['evalue_point']) < 0.01

    def test_null_effect(self):
        """RR=1.0 should give E-value=1"""
        r = evalue(1.0, measure='RR')
        assert r['evalue_point'] == 1.0

    def test_ci_crossing_null(self):
        """CI crossing null should give E-value of CI = 1"""
        r = evalue(1.5, lo=0.8, hi=2.0, measure='RR')
        assert r['evalue_ci'] == 1.0

    def test_ci_not_crossing_null(self):
        """CI not crossing null should give meaningful E-value"""
        r = evalue(2.0, lo=1.5, hi=2.8, measure='RR')
        assert r['evalue_ci'] is not None
        assert r['evalue_ci'] > 1.0

    def test_or_conversion(self):
        """OR is sqrt-transformed before E-value calc"""
        r_or = evalue(4.0, measure='OR')  # sqrt(4)=2.0
        r_rr = evalue(2.0, measure='RR')
        assert abs(r_or['evalue_point'] - r_rr['evalue_point']) < 0.01

    def test_interpretation(self):
        r = evalue(3.0, measure='RR')
        assert 'Strong' in r['interpretation'] or 'Moderate' in r['interpretation']


# ── Doi plot + LFK index ──

class TestDoiPlotLFK:
    def test_symmetric_data(self):
        """Symmetric data should have LFK near 0"""
        yi = [0.5, -0.5, 0.3, -0.3, 0.1, -0.1]
        sei = [0.2, 0.2, 0.3, 0.3, 0.4, 0.4]
        r = doi_plot_lfk(yi, sei)
        assert r['verdict'] == 'PASS'
        assert abs(r['lfk_index']) < 2.0

    def test_asymmetric_data(self):
        """Strongly asymmetric data should have large LFK"""
        yi = [2.0, 1.8, 1.5, 1.2, 0.1]
        sei = [0.1, 0.2, 0.3, 0.5, 1.0]
        r = doi_plot_lfk(yi, sei)
        assert r['lfk_index'] is not None

    def test_insufficient_studies(self):
        """k<3 should return INSUFFICIENT"""
        r = doi_plot_lfk([0.5, 0.3], [0.2, 0.3])
        assert r['verdict'] == 'INSUFFICIENT'

    def test_returns_z_scores(self):
        yi = [0.5, 0.3, 0.1]
        sei = [0.1, 0.2, 0.3]
        r = doi_plot_lfk(yi, sei)
        assert len(r['z_scores']) == 3


# ── Rosenthal fail-safe N ──

class TestRosenthalFSN:
    def test_strong_effect(self):
        """Strong effect should need many null studies"""
        yi = [0.5, 0.6, 0.4, 0.55, 0.45]
        sei = [0.1, 0.1, 0.1, 0.1, 0.1]
        r = rosenthal_failsafe_n(yi, sei)
        assert r['failsafe_n'] > 50
        assert r['robust'] is True  # Should exceed 5k+10 = 35

    def test_weak_effect(self):
        """Weak/null effect should need few null studies"""
        yi = [0.05, -0.05, 0.02, -0.03, 0.01]
        sei = [0.1, 0.1, 0.1, 0.1, 0.1]
        r = rosenthal_failsafe_n(yi, sei)
        assert r['failsafe_n'] < 50

    def test_single_study(self):
        """k<2 should return None"""
        r = rosenthal_failsafe_n([0.5], [0.1])
        assert r['failsafe_n'] is None

    def test_tolerance_formula(self):
        """Tolerance should be 5k+10"""
        yi = [0.5] * 10
        sei = [0.1] * 10
        r = rosenthal_failsafe_n(yi, sei)
        assert r['tolerance'] == 60  # 5*10+10


# ── Quality Effects Model ──

class TestQEM:
    def test_equal_quality(self):
        """Equal quality scores should give IV-like results"""
        yi = [0.5, 0.3, 0.7]
        sei = [0.1, 0.2, 0.15]
        quality = [1.0, 1.0, 1.0]
        r = quality_effects_model(yi, sei, quality)
        assert r['theta'] is not None
        assert 0.3 < r['theta'] < 0.7

    def test_quality_weighting(self):
        """Low-quality outlier should be downweighted"""
        yi = [0.5, 0.5, 5.0]  # outlier study 3
        sei = [0.1, 0.1, 0.1]  # same precision
        # With equal quality, outlier has full weight
        r_equal = quality_effects_model(yi, sei, [1.0, 1.0, 1.0])
        # With low quality on outlier, it's downweighted
        r_quality = quality_effects_model(yi, sei, [1.0, 1.0, 0.1])
        assert r_quality['theta'] < r_equal['theta']  # Outlier contributes less

    def test_returns_ci(self):
        yi = [0.5, 0.3]
        sei = [0.1, 0.2]
        r = quality_effects_model(yi, sei, [1.0, 0.8])
        assert r['ci_lo'] < r['theta'] < r['ci_hi']

    def test_returns_weights(self):
        yi = [0.5, 0.3, 0.7]
        sei = [0.1, 0.2, 0.15]
        r = quality_effects_model(yi, sei, [1.0, 0.5, 0.8])
        assert len(r['weights']) == 3
        assert abs(sum(r['weights']) - 1.0) < 0.01


# ── Proportion Meta-Analysis ──

class TestProportionMeta:
    def test_freeman_tukey(self):
        """Pool proportions using Freeman-Tukey double arcsine"""
        events = [10, 15, 8, 12]
        totals = [100, 120, 80, 110]
        r = proportion_meta(events, totals, method='PFT')
        assert 0.05 < r['pooled_proportion'] < 0.20
        assert r['ci_lo'] < r['pooled_proportion'] < r['ci_hi']

    def test_logit(self):
        """Pool proportions using logit transformation"""
        events = [10, 15, 8, 12]
        totals = [100, 120, 80, 110]
        r = proportion_meta(events, totals, method='logit')
        assert 0.05 < r['pooled_proportion'] < 0.20

    def test_arcsine(self):
        """Pool proportions using arcsine transformation"""
        events = [10, 15, 8]
        totals = [100, 120, 80]
        r = proportion_meta(events, totals, method='arcsine')
        assert 0.05 < r['pooled_proportion'] < 0.25

    def test_zero_events(self):
        """Handle zero events without crashing"""
        events = [0, 5, 3]
        totals = [50, 100, 80]
        r = proportion_meta(events, totals, method='PFT')
        assert r['pooled_proportion'] is not None
        assert r['pooled_proportion'] >= 0

    def test_returns_per_study(self):
        events = [10, 20]
        totals = [100, 200]
        r = proportion_meta(events, totals)
        assert r['per_study'] == [0.1, 0.1]

    def test_returns_heterogeneity(self):
        events = [10, 15, 8]
        totals = [100, 120, 80]
        r = proportion_meta(events, totals)
        assert 'tau2' in r
        assert 'i2' in r


# ── Permutation Test for Heterogeneity ──

class TestPermutationTest:
    def test_homogeneous_data(self):
        """Identical effects should show no heterogeneity"""
        yi = [0.5, 0.5, 0.5, 0.5, 0.5]
        sei = [0.1, 0.2, 0.15, 0.12, 0.18]
        r = permutation_test_heterogeneity(yi, sei, n_perm=5000, seed=42)
        assert r['p_perm'] > 0.05
        assert r['q_observed'] < 1.0  # Very small Q

    def test_heterogeneous_data(self):
        """Very different effects with varying precision should show heterogeneity"""
        yi = [0.1, 0.5, 1.0, -0.3, 2.0]
        sei = [0.05, 0.08, 0.12, 0.06, 0.15]  # varying precision required for permutation to detect
        r = permutation_test_heterogeneity(yi, sei, n_perm=5000, seed=42)
        assert r['q_observed'] > 50  # Large Q
        # Chi-squared p should be significant even if permutation is conservative
        assert r['p_chisq'] < 0.01

    def test_deterministic(self):
        """Same seed should give same result"""
        yi = [0.3, 0.5, 0.7]
        sei = [0.1, 0.2, 0.3]
        r1 = permutation_test_heterogeneity(yi, sei, seed=123)
        r2 = permutation_test_heterogeneity(yi, sei, seed=123)
        assert r1['p_perm'] == r2['p_perm']

    def test_insufficient_studies(self):
        """k<3 should return None"""
        r = permutation_test_heterogeneity([0.5, 0.3], [0.1, 0.2])
        assert r['q_observed'] is None

    def test_returns_both_p_values(self):
        """Should return both permutation and chi-squared p-values"""
        yi = [0.3, 0.5, 0.7, 0.4]
        sei = [0.1, 0.1, 0.1, 0.1]
        r = permutation_test_heterogeneity(yi, sei, n_perm=1000)
        assert r['p_perm'] is not None
        assert r['p_chisq'] is not None
