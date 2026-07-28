#!/usr/bin/env python3
"""
NACA 0012 CFD pre-analysis calculator.

Reproduces every hand calculation quoted in ../naca0012-airfoil/README.md from the
raw case inputs, so the numbers in the write-up can be checked rather than trusted:

  * chord Reynolds number and flow regime
  * free-stream dynamic pressure and inlet velocity components
  * thin-aerofoil lift estimate
  * sectional forces per unit span
  * y+ / first-cell-height sizing for either wall-function or wall-resolved meshing
  * verification helpers: Richardson extrapolation, GCI, mass-imbalance

Standard library only. Run with no arguments for the baseline case:

    python tools/preanalysis.py

Override any input from the command line, e.g. a wall-resolved mesh at 5 degrees:

    python tools/preanalysis.py --alpha 5 --y-plus 1

Jad El Badaoui - Aerospace Engineering, University of Bristol
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

# --------------------------------------------------------------------------------------
# Baseline case (Cornell MAE 5230 / ANSYS Fluent NACA 0012 module)
# --------------------------------------------------------------------------------------

CHORD = 1.00  # m
V_INF = 51.45  # m/s
ALPHA_DEG = 10.0  # degrees
RHO = 1.1767  # kg/m^3
MU = 1.009e-5  # kg/(m s)

# Typical sectional drag for this aerofoil, used only to illustrate the magnitude of the
# sectional forces in the pre-analysis. Reference: Ladson, NASA NACA 0012 validation data.
CD_REFERENCE = 0.012

# Log-law constants.
KAPPA, B_LOG = 0.41, 5.2


@dataclass
class Case:
    chord: float
    v_inf: float
    alpha_deg: float
    rho: float
    mu: float

    @property
    def alpha_rad(self) -> float:
        return math.radians(self.alpha_deg)

    @property
    def nu(self) -> float:
        """Kinematic viscosity, m^2/s."""
        return self.mu / self.rho

    @property
    def reynolds(self) -> float:
        """Chord Reynolds number: inertial transport / molecular momentum diffusion."""
        return self.rho * self.v_inf * self.chord / self.mu

    @property
    def q_inf(self) -> float:
        """Free-stream dynamic pressure, Pa."""
        return 0.5 * self.rho * self.v_inf**2

    @property
    def inlet_components(self) -> tuple[float, float]:
        """Velocity-inlet components for an aerofoil left aligned with the x-axis."""
        return (
            self.v_inf * math.cos(self.alpha_rad),
            self.v_inf * math.sin(self.alpha_rad),
        )

    @property
    def cl_thin_aerofoil(self) -> float:
        """Thin-aerofoil theory: CL = 2*pi*(alpha - alpha_0), with alpha_0 = 0 for a
        symmetric section. Inviscid and attached-flow - predicts no viscous drag."""
        return 2.0 * math.pi * self.alpha_rad

    def force_per_span(self, coefficient: float) -> float:
        """Sectional force per unit span, N/m."""
        return coefficient * self.q_inf * self.chord


# --------------------------------------------------------------------------------------
# Near-wall sizing
# --------------------------------------------------------------------------------------


def skin_friction_schlichting(re_x: float) -> float:
    """Local skin-friction coefficient for a turbulent flat-plate boundary layer.

    Schlichting correlation, valid for roughly 1e5 < Re_x < 1e9. Used only to get a
    first-cell height into the right order of magnitude before the first solve; the
    real y+ distribution must be plotted from the solution and the mesh then revised.
    """
    return 0.0576 * re_x ** (-0.2)


def first_cell_height(case: Case, y_plus_target: float, re_x: float | None = None) -> dict:
    """Estimate the wall-normal distance to the first cell centre for a target y+.

    Inverts  y+ = rho * y * u_tau / mu  after estimating wall shear from a flat-plate
    correlation:  tau_w = 0.5 * rho * V^2 * Cf,  u_tau = sqrt(tau_w / rho).
    """
    re_x = case.reynolds if re_x is None else re_x
    cf = skin_friction_schlichting(re_x)
    tau_w = 0.5 * case.rho * case.v_inf**2 * cf
    u_tau = math.sqrt(tau_w / case.rho)
    y1 = y_plus_target * case.mu / (case.rho * u_tau)
    return {"cf": cf, "tau_w": tau_w, "u_tau": u_tau, "y1": y1}


def inflation_total_thickness(y1: float, layers: int, growth: float) -> float:
    """Total height of an inflation stack: geometric series of `layers` terms."""
    if math.isclose(growth, 1.0):
        return y1 * layers
    return y1 * (growth**layers - 1.0) / (growth - 1.0)


def wall_law_u_plus(y_plus: float) -> tuple[float, float]:
    """Viscous-sublayer and log-layer velocity profiles at a given y+."""
    return y_plus, (1.0 / KAPPA) * math.log(y_plus) + B_LOG


def classify_y_plus(y_plus: float) -> str:
    if y_plus < 5:
        return "viscous sublayer - suits enhanced wall treatment / wall-resolved meshing"
    if y_plus < 30:
        return "BUFFER LAYER - avoid: neither the viscous nor the log-law fit applies"
    if y_plus <= 300:
        return "log layer - suits standard wall functions"
    return "beyond the log layer - first cell too far from the wall"


# --------------------------------------------------------------------------------------
# Verification helpers
# --------------------------------------------------------------------------------------


def grid_convergence_index(phi_1: float, phi_2: float, r: float, p: float = 2.0) -> dict:
    """Richardson extrapolation and GCI for two systematically refined meshes.

    phi_1 is the fine-mesh value, phi_2 the coarse-mesh value, r the refinement ratio
    and p the observed order of accuracy. Valid only for monotonic convergence with a
    consistent refinement ratio.
    """
    denominator = r**p - 1.0
    phi_ext = phi_1 + (phi_1 - phi_2) / denominator
    gci = 1.25 * abs((phi_1 - phi_2) / phi_1) / denominator * 100.0
    return {"phi_extrapolated": phi_ext, "gci_percent": gci}


def mass_imbalance(m_in: float, m_out: float) -> float:
    """Normalized mass-conservation error, %."""
    return abs(m_in - m_out) / m_in * 100.0


# --------------------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------------------


def rule(title: str) -> None:
    print(f"\n{title}\n{'-' * len(title)}")


def report(case: Case, y_plus_target: float, layers: int, growth: float) -> None:
    ux, uy = case.inlet_components

    print("=" * 74)
    print("  NACA 0012 - CFD PRE-ANALYSIS AND VALIDATION SUMMARY")
    print("=" * 74)

    rule("1. Case definition")
    print(f"  Chord, c                    {case.chord:>12.3f}  m")
    print(f"  Free-stream speed, V_inf    {case.v_inf:>12.2f}  m/s")
    print(f"  Angle of attack, alpha      {case.alpha_deg:>12.1f}  deg")
    print(f"  Density, rho                {case.rho:>12.4f}  kg/m^3")
    print(f"  Dynamic viscosity, mu       {case.mu:>12.3e}  kg/(m s)")
    print(f"  Kinematic viscosity, nu     {case.nu:>12.3e}  m^2/s")

    rule("2. Reynolds number and flow regime")
    print(f"  Re_c = rho*V*c/mu           {case.reynolds:>12.3e}")
    regime = "turbulent boundary layer and wake expected" if case.reynolds > 5e5 else "transitional or laminar"
    print(f"  Regime                      {regime}")
    print("  -> RANS with a turbulence closure is required; DNS is impractical here.")

    rule("3. Free-stream and inlet conditions")
    print(f"  Dynamic pressure, q_inf     {case.q_inf:>12.1f}  Pa")
    print(f"  Inlet U_x = V*cos(alpha)    {ux:>12.3f}  m/s")
    print(f"  Inlet U_y = V*sin(alpha)    {uy:>12.3f}  m/s")
    print(f"  Drag direction e_D          ({math.cos(case.alpha_rad):.4f}, {math.sin(case.alpha_rad):.4f})")
    print(f"  Lift direction e_L          ({-math.sin(case.alpha_rad):.4f}, {math.cos(case.alpha_rad):.4f})")

    rule("4. Thin-aerofoil hand prediction")
    cl_thin = case.cl_thin_aerofoil
    print(f"  CL = 2*pi*alpha             {cl_thin:>12.3f}")
    print(f"  Lift per span, L'           {case.force_per_span(cl_thin):>12.1f}  N/m")
    print(f"  Drag per span at CD={CD_REFERENCE:<5g}  {case.force_per_span(CD_REFERENCE):>12.1f}  N/m")
    print("  -> Inviscid theory predicts zero drag, so CFD and experiment remain necessary.")

    rule(f"5. Near-wall sizing for y+ = {y_plus_target:g}")
    nw = first_cell_height(case, y_plus_target)
    total = inflation_total_thickness(nw["y1"], layers, growth)
    u_visc, u_log = wall_law_u_plus(y_plus_target)
    print(f"  Skin friction, C_f          {nw['cf']:>12.5f}   (flat-plate estimate)")
    print(f"  Wall shear stress, tau_w    {nw['tau_w']:>12.3f}  Pa")
    print(f"  Friction velocity, u_tau    {nw['u_tau']:>12.4f}  m/s")
    print(f"  First cell height, y_1      {nw['y1']:>12.3e}  m   ({nw['y1'] / case.chord:.3e} c)")
    print(f"  Inflation stack             {layers} layers, growth {growth:g}")
    print(f"  Total inflation thickness   {total:>12.3e}  m   ({total / case.chord:.4f} c)")
    print(f"  u+ viscous law / log law    {u_visc:>12.2f} / {u_log:.2f}")
    print(f"  Placement                   {classify_y_plus(y_plus_target)}")

    rule("6. Experimental reference conditions")
    print(f"  Reference source            NASA NACA 0012 validation data")
    print(f"    Surface Cp                Gregory and O'Reilly")
    print(f"    Force coefficients        Ladson")
    print(f"  Match before comparing      Re, angle of attack, reference area/chord definitions")
    print()
    print("  Compare the surface Cp distribution over the complete chord, not just an")
    print("  integrated coefficient: agreement in a single scalar can arise from error")
    print("  cancellation, whereas a point-by-point distribution match cannot.")

    rule("7. Worked verification examples")
    gci = grid_convergence_index(phi_1=2.4180, phi_2=2.3950, r=2.0, p=2.0)
    print(f"  GCI example, generic quantity (fine 2.4180, coarse 2.3950, r=2, p=2)")
    print(f"    Richardson-extrapolated value {gci['phi_extrapolated']:.5f}")
    print(f"    GCI_12                        {gci['gci_percent']:.2f}%")
    print(f"  Mass-imbalance check (1.0000000 in, 0.9999999 out)"
          f"   {mass_imbalance(1.0, 0.9999999):.1e}%")

    print("\n" + "=" * 74)
    print("  A converged solution is not a validated one. See naca0012-airfoil/README.md")
    print("=" * 74)


def main() -> None:
    p = argparse.ArgumentParser(
        description="NACA 0012 CFD pre-analysis calculator.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--chord", type=float, default=CHORD, help="chord length, m")
    p.add_argument("--v-inf", type=float, default=V_INF, help="free-stream speed, m/s")
    p.add_argument("--alpha", type=float, default=ALPHA_DEG, help="angle of attack, deg")
    p.add_argument("--rho", type=float, default=RHO, help="density, kg/m^3")
    p.add_argument("--mu", type=float, default=MU, help="dynamic viscosity, kg/(m s)")
    p.add_argument("--y-plus", type=float, default=30.0, help="target first-cell y+")
    p.add_argument("--layers", type=int, default=10, help="inflation layer count")
    p.add_argument("--growth", type=float, default=1.2, help="inflation growth rate")
    args = p.parse_args()

    case = Case(
        chord=args.chord,
        v_inf=args.v_inf,
        alpha_deg=args.alpha,
        rho=args.rho,
        mu=args.mu,
    )
    report(case, args.y_plus, args.layers, args.growth)


if __name__ == "__main__":
    main()
