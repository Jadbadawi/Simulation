"""Tests for the NACA 0012 pre-analysis calculator.

The point of `preanalysis.py` is that the numbers quoted in the write-ups can be
checked rather than trusted. These tests pin the quoted values, so if the calculator
changes, the mismatch with the README shows up here instead of going unnoticed.

Standard library only, like the module under test:

    python -m pytest tools/
"""

from __future__ import annotations

import math

import pytest

import preanalysis as pa


@pytest.fixture
def baseline() -> pa.Case:
    """The baseline case documented in naca0012-airfoil/README.md."""
    return pa.Case(
        chord=pa.CHORD,
        v_inf=pa.V_INF,
        alpha_deg=pa.ALPHA_DEG,
        rho=pa.RHO,
        mu=pa.MU,
    )


# --------------------------------------------------------------------------------------
# Values quoted in the write-up
# --------------------------------------------------------------------------------------

def test_reynolds_number_is_six_million(baseline):
    """The case is set up to hit Re_c = 6e6, matching the NASA validation data."""
    assert baseline.reynolds == pytest.approx(6.0e6, rel=1e-3)


def test_dynamic_pressure(baseline):
    """q = 0.5 * rho * V^2, quoted as 1557.4 Pa."""
    assert baseline.q_inf == pytest.approx(1557.4, rel=1e-3)
    assert baseline.q_inf == pytest.approx(0.5 * baseline.rho * baseline.v_inf**2)


def test_inlet_velocity_components(baseline):
    """Inlet decomposed at 10 degrees incidence: 50.668 and 8.934 m/s."""
    ux, uy = baseline.inlet_components
    assert ux == pytest.approx(50.668, abs=1e-3)
    assert uy == pytest.approx(8.934, abs=1e-3)


def test_inlet_components_preserve_magnitude(baseline):
    ux, uy = baseline.inlet_components
    assert math.hypot(ux, uy) == pytest.approx(baseline.v_inf)


def test_thin_aerofoil_lift_estimate(baseline):
    """2*pi*alpha, quoted as 1.097 -- computed before any solver was opened."""
    assert baseline.cl_thin_aerofoil == pytest.approx(1.097, abs=1e-3)
    assert baseline.cl_thin_aerofoil == pytest.approx(2 * math.pi * baseline.alpha_rad)


def test_first_cell_height_for_wall_functions(baseline):
    """y+ = 30 sizing: u_tau = 1.8334 m/s and y1 = 1.403e-4 m."""
    result = pa.first_cell_height(baseline, y_plus_target=30.0)
    assert result["u_tau"] == pytest.approx(1.8334, abs=1e-3)
    assert result["y1"] == pytest.approx(1.403e-4, rel=1e-2)


# --------------------------------------------------------------------------------------
# Physical relationships
# --------------------------------------------------------------------------------------

def test_kinematic_viscosity(baseline):
    assert baseline.nu == pytest.approx(baseline.mu / baseline.rho)


def test_zero_incidence_gives_zero_lift(baseline):
    """A symmetric aerofoil at zero incidence produces no lift."""
    case = pa.Case(baseline.chord, baseline.v_inf, 0.0, baseline.rho, baseline.mu)
    assert case.cl_thin_aerofoil == pytest.approx(0.0)
    ux, uy = case.inlet_components
    assert ux == pytest.approx(baseline.v_inf)
    assert uy == pytest.approx(0.0)


def test_first_cell_height_scales_inversely_with_y_plus(baseline):
    """Halving the y+ target should halve the first cell height."""
    coarse = pa.first_cell_height(baseline, y_plus_target=30.0)["y1"]
    fine = pa.first_cell_height(baseline, y_plus_target=15.0)["y1"]
    assert fine == pytest.approx(coarse / 2.0, rel=1e-9)


def test_wall_resolved_mesh_is_much_finer(baseline):
    """y+ = 1 demands a far finer first cell than a wall-function mesh."""
    wall_fn = pa.first_cell_height(baseline, y_plus_target=30.0)["y1"]
    resolved = pa.first_cell_height(baseline, y_plus_target=1.0)["y1"]
    assert resolved < wall_fn / 20.0


def test_skin_friction_decreases_with_reynolds():
    """Cf falls as Re rises for a turbulent boundary layer."""
    assert pa.skin_friction_schlichting(1e6) > pa.skin_friction_schlichting(1e7)


def test_inflation_total_thickness_matches_geometric_series():
    y1, layers, growth = 1e-4, 10, 1.2
    expected = y1 * (growth**layers - 1) / (growth - 1)
    assert pa.inflation_total_thickness(y1, layers, growth) == pytest.approx(expected)


def test_inflation_thickness_grows_with_layer_count():
    a = pa.inflation_total_thickness(1e-4, 5, 1.2)
    b = pa.inflation_total_thickness(1e-4, 15, 1.2)
    assert b > a


# --------------------------------------------------------------------------------------
# Near-wall classification
# --------------------------------------------------------------------------------------

def test_wall_law_in_viscous_sublayer():
    """In the sublayer u+ = y+."""
    linear, _ = pa.wall_law_u_plus(2.0)
    assert linear == pytest.approx(2.0)


def test_wall_law_log_layer_matches_log_law():
    _, log_value = pa.wall_law_u_plus(100.0)
    assert log_value == pytest.approx(math.log(100.0) / pa.KAPPA + pa.B_LOG)


@pytest.mark.parametrize("y_plus", [0.5, 1.0, 5.0, 15.0, 30.0, 100.0, 500.0])
def test_classify_y_plus_returns_a_description(y_plus):
    assert isinstance(pa.classify_y_plus(y_plus), str)
    assert pa.classify_y_plus(y_plus).strip() != ""


def test_classify_y_plus_distinguishes_regimes():
    """A wall-resolved and a wall-function placement must not be described alike."""
    assert pa.classify_y_plus(1.0) != pa.classify_y_plus(30.0)


# --------------------------------------------------------------------------------------
# Verification helpers
# --------------------------------------------------------------------------------------

def test_grid_convergence_index_zero_for_identical_solutions():
    """No change between grids means no discretisation error to report."""
    result = pa.grid_convergence_index(1.0, 1.0, r=2.0)
    assert result["gci_percent"] == pytest.approx(0.0)


def test_richardson_extrapolation_overshoots_finer_grid():
    """With monotone convergence the extrapolated value lies beyond the fine grid."""
    phi_1, phi_2 = 1.02, 1.05  # fine, coarse
    result = pa.grid_convergence_index(phi_1, phi_2, r=2.0, p=2.0)
    assert result["phi_extrapolated"] < phi_1
    assert result["gci_percent"] > 0.0


def test_gci_shrinks_as_grids_converge():
    """A smaller gap between grids implies a smaller error band."""
    far = pa.grid_convergence_index(1.00, 1.20, r=2.0)["gci_percent"]
    near = pa.grid_convergence_index(1.00, 1.01, r=2.0)["gci_percent"]
    assert near < far


def test_mass_imbalance_is_zero_when_conserved():
    assert pa.mass_imbalance(10.0, 10.0) == pytest.approx(0.0)


def test_mass_imbalance_is_normalised_by_inflow():
    """Reported as a percentage, so a 1% shortfall reads 1.0 at any flow rate."""
    assert pa.mass_imbalance(100.0, 99.0) == pytest.approx(1.0)
    assert pa.mass_imbalance(1.0, 0.99) == pytest.approx(1.0)


# --------------------------------------------------------------------------------------
# The documented command still works
# --------------------------------------------------------------------------------------

def test_report_runs_without_error(baseline, capsys):
    pa.report(baseline, y_plus_target=30.0, layers=10, growth=1.2)
    out = capsys.readouterr().out
    assert "Reynolds number" in out
    assert "6.000e+06" in out
