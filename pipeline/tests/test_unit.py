"""Pruebas unitarias de las piezas del pipeline."""

import numpy as np
import pytest
from scipy import ndimage

from instrastro import alignment, detection, photometry, reduction, seeing
from instrastro.synthetic import _render


def _star_field(size=400, n=60, fwhm=6.0, seed=1, sky=20.0):
    rng = np.random.default_rng(seed)
    img = np.full((size, size), sky)
    xs, ys = rng.uniform(30, size - 30, n), rng.uniform(30, size - 30, n)
    for x, y in zip(xs, ys):
        _render(img, x, y, rng.uniform(2e3, 2e4), fwhm / 2.3548)
    return img + rng.normal(0, 1.0, img.shape), xs, ys


def test_tycho_to_johnson_hip116375():
    # HIP 116375: BT = 8.453, VT = 7.191 -> V = 7.08, B-V = 1.07 (SIMBAD / Hipparcos)
    b, v = photometry.tycho_to_johnson(8.453, 7.191)
    assert v == pytest.approx(7.077, abs=0.002)
    assert b - v == pytest.approx(1.073, abs=0.002)


def test_combine_rejects_cosmic_ray():
    rng = np.random.default_rng(0)
    frames = [rng.normal(100, 1, (50, 50)).astype(np.float32) for _ in range(7)]
    frames[3][10, 10] = 5000.0
    frames[2][20, 20] = np.nan
    out = reduction.combine(frames, method="mean")
    assert out[10, 10] == pytest.approx(100, abs=3)
    assert np.isfinite(out[20, 20])


def test_measure_shift_subpixel():
    img, _, _ = _star_field()
    shifted = ndimage.shift(img, (3.3, -5.6), order=3)
    dy, dx = alignment.measure_shift(img, shifted)
    assert dy == pytest.approx(-3.3, abs=0.1)
    assert dx == pytest.approx(5.6, abs=0.1)
    realigned = alignment.shift_image(shifted.astype(np.float32), dy, dx)
    core = (slice(20, -20), slice(20, -20))
    assert np.nanmax(np.abs(realigned[core] - img[core])) < 0.05 * img.max()


def test_gaussian_fit_recovers_fwhm():
    img = np.zeros((41, 41))
    _render(img, 20.3, 19.6, 1e4, 5.0 / 2.3548)
    fit = seeing.fit_gaussian(img, 2.0)
    assert fit["fwhm"] == pytest.approx(5.0, rel=0.02)
    assert fit["x0"] == pytest.approx(20.3, abs=0.05)


def test_measure_psf_on_field():
    img, _, _ = _star_field(fwhm=6.0)
    psf = seeing.measure_psf(img.astype(np.float32), fwhm_guess_px=8.0, n_sigma=5)
    summary = seeing.summarize(psf, 0.36)
    assert summary["n_estrellas"] >= 5
    assert summary["fwhm_px"] == pytest.approx(6.0, rel=0.05)


def test_diffraction_limit_mas500():
    d = seeing.diffraction_limit("V", 0.5)
    assert d["rayleigh_arcsec"] == pytest.approx(0.277, abs=0.005)


def test_aperture_photometry_and_correction():
    rng = np.random.default_rng(3)
    fwhm = 7.0
    img = np.full((512, 512), 15.0)
    grid = np.arange(60, 500, 90.0)
    x = np.repeat(grid, len(grid)) + rng.uniform(-0.5, 0.5, len(grid) ** 2)
    y = np.tile(grid, len(grid)) + rng.uniform(-0.5, 0.5, len(grid) ** 2)
    for xi, yi in zip(x, y):
        _render(img, xi, yi, 2000.0, fwhm / 2.3548)
    img = rng.poisson(img * 150) / 150.0
    flux, err, _ = photometry.aperture_photometry_list(img, x, y, fwhm, 3.5 * fwhm, 5 * fwhm, gain_eff=150)
    # Una gaussiana tiene el 93.75 % del flujo dentro de r = FWHM
    assert np.median(flux) == pytest.approx(2000 * 0.9375, rel=0.01)
    corr, corr_err, n = photometry.aperture_correction(img, x, y, fwhm, 3 * fwhm, 3.5 * fwhm, 5 * fwhm)
    assert corr == pytest.approx(1 / 0.9375, rel=0.01)
    assert n >= 3


def test_zeropoint_roundtrip():
    zp, k, x_air, m = 21.2, 0.15, 1.3, 7.08
    flux = 10 ** (-0.4 * (m - zp + k * x_air))
    assert photometry.zeropoint(flux, m, k, x_air) == pytest.approx(zp)


def test_crossmatch_one_to_one():
    from astropy.table import Table
    tb = Table({"x": [10.0, 10.5, 50.0], "y": [10.0, 10.2, 50.0], "mag": [15.0, 15.5, 16.0],
                "mag_err": [0.01] * 3, "saturated": [False] * 3})
    tv = Table({"x": [10.2, 50.1], "y": [10.1, 49.9], "mag": [14.0, 15.2], "mag_err": [0.01] * 2,
                "saturated": [False] * 2})
    out = detection.crossmatch(tb, tv, tol_px=2.0)
    assert len(out) == 2
    assert sorted(np.round(out["B_V"], 2)) == [0.8, 1.0]


def test_limiting_magnitude_monotonic():
    assert detection.limiting_magnitude(0.5, 8, 21) < detection.limiting_magnitude(0.25, 8, 21)


def test_color_term_fit_and_correction():
    """Con c = 2,3 el ZP depende del color; el catálogo B corregido debe recuperar B y B−V."""
    from astropy.table import Table
    zp0, c = 19.43, 2.3
    rows = []
    for name, bv in (("azul", 0.05), ("roja", 0.28)):
        for i in range(5):
            rows.append((name, "B", 1.0, False, bv, zp0 + c * bv + 0.002 * (i - 2)))
            rows.append((name, "V", 1.0, False, bv, 21.75))
    tbl = Table(rows=rows, names=("estrella", "filtro", "airmass", "saturada", "B_V", "zp"))
    res = photometry.summarize_zeropoints(tbl, color_terms=["B"])
    b, v = res[res["filtro"] == "B"][0], res[res["filtro"] == "V"][0]
    assert b["zp"] == pytest.approx(zp0, abs=1e-3) and b["c"] == pytest.approx(c, abs=1e-2)
    assert 0 < b["c_err"] < 0.2
    assert v["c"] == 0 and v["zp"] == pytest.approx(21.75)
    # Sin término de color, B queda como un solo ZP promedio
    assert photometry.summarize_zeropoints(tbl)["c"][0] == 0

    # Estrella con B−V = 0,8: el catálogo B calibrado con ZP0 da B0 = B − c (B−V)
    true_b, true_v = 15.8, 15.0
    cmd = Table({"B": [true_b - c * 0.8], "B_err": [0.02], "V": [true_v], "V_err": [0.02]})
    detection.apply_color_term(cmd, c, 0.1)
    assert cmd["B_V"][0] == pytest.approx(0.8) and cmd["B"][0] == pytest.approx(true_b)
    assert cmd["B0"][0] == pytest.approx(true_b - c * 0.8)
