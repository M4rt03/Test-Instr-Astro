"""Prueba de punta a punta: el pipeline debe recuperar los valores verdaderos de un set sintético."""

import json

import numpy as np
import pytest
from astropy.table import Table
from scipy.spatial import cKDTree

from instrastro import pipeline, synthetic


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    outdir = tmp_path_factory.mktemp("demo")
    truth, config = synthetic.make_dataset(outdir)
    results = pipeline.run(config)
    return truth, results


def _match(tbl, truth, tol=3.0):
    stars = truth["stars"]
    dist, idx = cKDTree(np.column_stack([stars["x"], stars["y"]])).query(
        np.column_stack([tbl["x"], tbl["y"]]), distance_upper_bound=tol)
    ok = np.isfinite(dist)
    return ok, np.where(ok, idx, 0)


def test_calibration_frames(demo):
    truth, results = demo
    cal = json.loads((results / "noche1" / "calibracion.json").read_text())
    assert cal["dark_current_mediana_adu_s"] == pytest.approx(truth["dark_current"], rel=0.05)
    assert cal["ruido_lectura_adu"] == pytest.approx(truth["readnoise"], rel=0.1)


def test_zeropoints(demo):
    truth, results = demo
    zps = Table.read(results / "noche1" / "zeropoints.ecsv")
    for row in zps:
        # El flat se normaliza a la mediana central (~0.98 del centro óptico), lo que cambia el ZP
        # en ~0.02 mag; ese factor es común a estándares y ciencia y se cancela en las magnitudes.
        assert row["zp"] == pytest.approx(truth["zp"][row["filtro"]], abs=0.05)
        assert row["n"] == 6


def test_seeing(demo):
    truth, results = demo
    for night, data in truth["nights"].items():
        s = json.loads((results / night / "seeing.json").read_text())
        expected = data["fwhm_px"] * truth["pixscale"]
        for filt in ("B", "V"):
            assert s[filt]["seeing_exposiciones_mediana"] == pytest.approx(expected, rel=0.05)
            # El apilado se ensancha un poco por la interpolación del alineamiento
            assert s[filt]["seeing_apilado_arcsec"] == pytest.approx(expected, rel=0.1)


def test_catalog_magnitudes(demo):
    truth, results = demo
    for night in truth["nights"]:
        for filt in ("B", "V"):
            tbl = Table.read(results / night / f"catalogo_agresiva_{filt}.ecsv")
            ok, idx = _match(tbl, truth)
            ok &= ~np.asarray(tbl["saturated"])
            diff = np.asarray(tbl["mag"])[ok] - truth["stars"][filt][idx[ok]]
            err = np.asarray(tbl["mag_err"])[ok]
            assert ok.sum() >= 25
            assert np.median(diff) == pytest.approx(0.0, abs=0.03), f"{night} {filt}"
            # Los errores reportados deben ser realistas: residuos normalizados con dispersión <~ 1
            z = diff / err
            assert 1.4826 * np.median(np.abs(z - np.median(z))) < 1.5, f"{night} {filt}"
            bright = truth["stars"][filt][idx[ok]] < 15
            assert 1.4826 * np.median(np.abs(diff[bright] - np.median(diff[bright]))) < 0.06, f"{night} {filt}"


def test_color_index(demo):
    truth, results = demo
    tbl = Table.read(results / "noche1" / "catalogo_BV.ecsv")
    ok, idx = _match(tbl, truth)
    ok &= ~np.asarray(tbl["saturated"])
    diff = np.asarray(tbl["B_V"])[ok] - truth["stars"]["BV"][idx[ok]]
    assert np.median(diff) == pytest.approx(0.0, abs=0.04)


def test_saturated_star_flagged(demo):
    truth, results = demo
    tbl = Table.read(results / "noche1" / "catalogo_agresiva_V.ecsv")
    x0, y0 = truth["stars"]["x"][0], truth["stars"]["y"][0]
    near = np.hypot(tbl["x"] - x0, tbl["y"] - y0) < 5
    assert near.any() and np.all(tbl["saturated"][near])


def test_night_comparison(demo):
    _, results = demo
    tbl = Table.read(results / "comparacion_noches.ecsv")
    v = {row["noche"]: row for row in tbl if row["filtro"] == "V"}
    # La noche 2 (Luna, peor seeing) debe tener cielo más brillante y magnitud límite menor
    assert v["noche2"]["cielo_mag_arcsec2"] < v["noche1"]["cielo_mag_arcsec2"]
    assert v["noche2"]["mag_limite_5sigma"] < v["noche1"]["mag_limite_5sigma"]
    assert v["noche2"]["seeing_apilado"] > v["noche1"]["seeing_apilado"]


def test_outputs_exist(demo):
    _, results = demo
    for name in ("calibracion.png", "reduccion_B.png", "reduccion_V.png", "seeing.png", "segmentacion.png",
                 "color_magnitud.png", "color_rgb.png", "mapa_color_BV.png", "visibilidad.png"):
        assert (results / "noche1" / "figuras" / name).exists(), name
    for name in ("stack_B.fits", "stack_V.fits", "bitacora.tex", "estandares.tex", "exposiciones.tex"):
        assert (results / "noche1" / name).exists(), name
