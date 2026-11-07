#!/usr/bin/env python3
"""Validate the supplied core-shell case against an independent Mie solution.

Run after producing a fresh coreshell.SIFlux.AVSCAT file, for example
    python3 ValidateCoreShell.py CoreShellMie.SIFlux.AVSCAT

The case has core/outer radii 0.5/1 um, core epsilon 4+1j, shell epsilon
10+0.5j, vacuum host, and mu=1 throughout. Mie coefficients are computed
by matching the Riccati-Bessel fields at both interfaces. Core absorption
is the inward Poynting flux through the inner interface. Shell absorption
is the difference between the outer and inner inward fluxes.

Requires NumPy and SciPy. The 5% default is a BEM discretization tolerance,
not the tolerance for comparing a binary with its pre-cleanup output.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import spherical_jn, spherical_yn


def riccati(order, z, kind):
    function = spherical_jn if kind == "j" else spherical_yn
    value = function(order, z)
    return z * value, value + z * function(order, z, derivative=True)


def coated_sphere(omega, core_epsilon=4 + 1j, shell_epsilon=10 + 0.5j,
                  core_radius=0.5, outer_radius=1.0, order_limit=20):
    """Return cross sections in um^2 and the dimensionless asymmetry g."""
    core_index = np.sqrt(complex(core_epsilon))
    shell_index = np.sqrt(complex(shell_epsilon))
    coefficients = []
    core_flux = 0.0
    outer_flux = 0.0
    orders = np.arange(1, order_limit + 1)
    for order in orders:
        host, host_derivative = riccati(order, omega * outer_radius, "j")
        host_y, host_y_derivative = riccati(order, omega * outer_radius, "y")
        outgoing = host + 1j * host_y
        outgoing_derivative = host_derivative + 1j * host_y_derivative
        shell_outer, shell_outer_derivative = riccati(
            order, shell_index * omega * outer_radius, "j")
        shell_outer_y, shell_outer_y_derivative = riccati(
            order, shell_index * omega * outer_radius, "y")
        shell_inner, shell_inner_derivative = riccati(
            order, shell_index * omega * core_radius, "j")
        shell_inner_y, shell_inner_y_derivative = riccati(
            order, shell_index * omega * core_radius, "y")
        core, core_derivative = riccati(order, core_index * omega * core_radius, "j")
        pair = []
        for core_factor, shell_factor in ((1 / core_index, 1 / shell_index),
                                          (core_index, shell_index)):
            matrix = np.array([
                [outgoing, shell_outer, shell_outer_y, 0],
                [outgoing_derivative, shell_factor * shell_outer_derivative,
                 shell_factor * shell_outer_y_derivative, 0],
                [0, shell_inner, shell_inner_y, -core],
                [0, shell_factor * shell_inner_derivative,
                 shell_factor * shell_inner_y_derivative, -core_factor * core_derivative],
            ], dtype=complex)
            coefficient, _, _, core_amplitude = np.linalg.solve(
                matrix, [host, host_derivative, 0, 0])
            pair.append(coefficient)
            weight = 2 * order + 1
            inner_field = core_amplitude * core
            inner_derivative = core_factor * core_amplitude * core_derivative
            core_flux -= weight * np.imag(np.conj(inner_field) * inner_derivative)
            outer_field = host - coefficient * outgoing
            outer_derivative = host_derivative - coefficient * outgoing_derivative
            outer_flux -= weight * np.imag(np.conj(outer_field) * outer_derivative)
        coefficients.append(pair)

    a, b = np.array(coefficients).T
    prefactor = 2 * np.pi / omega**2
    scattering = prefactor * np.sum((2 * orders + 1) * (abs(a)**2 + abs(b)**2))
    extinction = prefactor * np.sum((2 * orders + 1) * np.real(a + b))
    absorption = extinction - scattering
    moment = np.sum((2 * orders + 1) / (orders * (orders + 1)) * np.real(a * np.conj(b)))
    moment += np.sum(orders[:-1] * (orders[:-1] + 2) / (orders[:-1] + 1)
                     * np.real(a[:-1] * np.conj(a[1:]) + b[:-1] * np.conj(b[1:])))
    moment *= 4 * np.pi / omega**2
    np.testing.assert_allclose(prefactor * outer_flux, absorption, rtol=1e-10, atol=1e-12)
    core_absorption = prefactor * core_flux
    return dict(Cabs=float(absorption), Csca=float(scattering),
                gCsca=float(moment), g=float(moment / scattering),
                Cabs_core=float(core_absorption), Cabs_shell=float(absorption - core_absorption))


def check_reference():
    """Check truncation, the homogeneous limit, and the quasistatic limit."""
    for omega in (0.3, 1.0):
        first = coated_sphere(omega, order_limit=12)
        second = coated_sphere(omega, order_limit=24)
        np.testing.assert_allclose(list(first.values()), list(second.values()), rtol=1e-11)
        index = np.sqrt(4 + 1j)
        orders = np.arange(1, 21)
        host, derivative = riccati(orders, omega, "j")
        host_y, derivative_y = riccati(orders, omega, "y")
        internal, internal_derivative = riccati(orders, index * omega, "j")
        outgoing, outgoing_derivative = host + 1j * host_y, derivative + 1j * derivative_y
        a = ((index * internal * derivative - host * internal_derivative)
             / (index * internal * outgoing_derivative - outgoing * internal_derivative))
        b = ((internal * derivative - index * host * internal_derivative)
             / (internal * outgoing_derivative - index * outgoing * internal_derivative))
        expected_sca = 2 * np.pi / omega**2 * np.sum(
            (2 * orders + 1) * (abs(a)**2 + abs(b)**2))
        expected_abs = 2 * np.pi / omega**2 * np.sum(
            (2 * orders + 1) * np.real(a + b)) - expected_sca
        homogeneous = coated_sphere(omega, shell_epsilon=4 + 1j)
        np.testing.assert_allclose([homogeneous["Cabs"], homogeneous["Csca"]],
                                   [expected_abs, expected_sca], rtol=1e-11)
        np.testing.assert_allclose(coated_sphere(omega, shell_epsilon=10)["Cabs_shell"],
                                   0, atol=1e-12)
        np.testing.assert_allclose(coated_sphere(omega, core_epsilon=4)["Cabs_core"],
                                   0, atol=1e-12)
    epsilon_core, epsilon_shell, fraction = 4 + 1j, 10 + 0.5j, 0.5**3
    polarizability = (((epsilon_shell - 1) * (epsilon_core + 2 * epsilon_shell)
                       + fraction * (epsilon_core - epsilon_shell) * (1 + 2 * epsilon_shell))
                      / ((epsilon_shell + 2) * (epsilon_core + 2 * epsilon_shell)
                         + 2 * fraction * (epsilon_shell - 1) * (epsilon_core - epsilon_shell)))
    omega = 1e-3
    small = coated_sphere(omega, order_limit=6)
    np.testing.assert_allclose(small["Cabs"], 4 * np.pi * omega * polarizability.imag, rtol=2e-5)
    np.testing.assert_allclose(small["Csca"], 8 * np.pi / 3 * omega**4 * abs(polarizability)**2,
                               rtol=2e-5)


def validate(path, tolerance):
    groups = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        row = line.split()
        if len(row) != 6:
            raise ValueError(f"{path}: expected six columns, got {row}")
        tag, omega, surface = row[:3]
        if tag != "0.0" or surface not in ("01", "02"):
            raise ValueError(f"{path}: unexpected transform or surface in {row}")
        values = np.array(row[3:], dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{path}: non-finite values in {row}")
        entries = groups.setdefault(float(omega), {})
        if surface in entries:
            raise ValueError(f"{path}: duplicate surface {surface} at omega {omega}")
        entries[surface] = values
    if set(groups) != {0.3, 1.0}:
        raise ValueError(f"{path}: expected frequencies 0.3 and 1.0, got {sorted(groups)}")
    passed = True
    report = []
    for omega, entries in sorted(groups.items()):
        if set(entries) != {"01", "02"}:
            raise ValueError(f"{path}: expected outer and inner surfaces at omega {omega}")
        total = entries["01"] + entries["02"]
        observed = dict(zip(("Cabs", "Csca", "gCsca"), map(float, total)))
        observed["Cabs_shell"] = float(entries["01"][0])
        observed["Cabs_core"] = float(entries["02"][0])
        expected = coated_sphere(omega)
        physical = (observed["Cabs_core"] >= 0 and observed["Cabs_shell"] >= 0
                    and total[1] > 0 and abs(total[2]) <= total[1])
        observed["g"] = float(total[2] / total[1]) if total[1] != 0 else float("nan")
        errors = {key: abs(value - expected[key]) / abs(expected[key])
                  for key, value in observed.items()}
        ok = physical and max(errors.values()) <= tolerance
        passed &= ok
        report.append(dict(omega=omega, passed=bool(ok), physical=bool(physical),
                           observed=observed, mie=expected, relative_errors=errors))
        print(f"{'PASS' if ok else 'FAIL'} {path.name} omega={omega:g} "
              f"max_relative_error={max(errors.values()):.3%} physical={physical}")
        for key in ("Cabs_shell", "Cabs_core", "Cabs", "Csca", "gCsca", "g"):
            print(f"  {key:12s} BEM={observed[key]:.10g} Mie={expected[key]:.10g}")
    return passed, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument("--rtol", type=float, default=0.05)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    if not 0 < args.rtol < 1:
        parser.error("--rtol must be between 0 and 1")
    check_reference()
    passed = True
    report = {}
    for path in args.files:
        ok, results = validate(path, args.rtol)
        passed &= ok
        report[str(path)] = results
    if args.json:
        args.json.write_text(json.dumps(report, indent=2) + "\n")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
