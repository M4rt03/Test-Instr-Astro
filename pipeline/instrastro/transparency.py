"""Corrección de transparencia (nubes) con estrellas del campo de un catálogo (APASS DR9).

Cuando el zeropoint viene de estándares de otra noche, una noche con nubes queda con todas las
magnitudes desplazadas en la misma cantidad (las nubes son grises). Se mide ese desplazamiento
comparando las magnitudes ya calibradas (incluido el término de color) de estrellas del campo con
las de APASS:

    Δ_V = mediana(V_APASS − V),     Δ_B = mediana(B_APASS − B)

y se corrigen los zeropoints de la noche. Con término de color en B (B = B0 + c (B−V)), para que
B y V finales se muevan en Δ_B y Δ_V hay que mover el ZP de B−V = 0 en

    δZP_B = Δ_B (1 − c) + c Δ_V          (= Δ si Δ_B = Δ_V = Δ)

La corrección también absorbe el error de la extinción supuesta y las diferencias de escala con
APASS (~0,02–0,05 mag); por eso es una corrección **de la noche**, no un zeropoint nuevo del sistema.
"""

import csv
import io
import logging
import urllib.parse
import urllib.request

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.table import Table

log = logging.getLogger("instrastro")

# Servidores de VizieR en orden de preferencia (CDS a veces no responde; el espejo chino solo por http)
SERVERS = ("https://vizier.cds.unistra.fr", "https://vizier.cfa.harvard.edu", "http://vizier.china-vo.org")
APASS = "II/336/apass9"


def query_apass(center, radius_arcmin, servers=SERVERS, timeout=60):
    """Estrellas de APASS DR9 con B y V dentro de ``radius_arcmin`` de ``center`` (SkyCoord).

    Usa la interfaz ASU de VizieR (texto separado por tabuladores) para no depender de astroquery.
    Prueba los servidores en orden y devuelve una tabla con ra, dec, B, V, eB, eV.
    """
    params = {"-source": APASS, "-c": f"{center.ra.deg:.6f} {center.dec.deg:+.6f}",
              "-c.rm": f"{radius_arcmin:.2f}", "-out": "RAJ2000,DEJ2000,Vmag,e_Vmag,Bmag,e_Bmag",
              "-out.max": "unlimited", "-oc.form": "d"}
    errors = []
    for server in servers:
        url = f"{server}/viz-bin/asu-tsv?{urllib.parse.urlencode(params)}"
        req = urllib.request.Request(url, headers={"User-Agent": "instrastro/1.0 (pipeline FIS2421)"})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                text = r.read().decode("utf-8", "replace")
        except Exception as e:      # noqa: BLE001 — se prueba el siguiente servidor
            errors.append(f"{server}: {e}")
            continue
        rows = [line for line in text.splitlines() if line and not line.startswith("#")]
        if len(rows) < 3:
            errors.append(f"{server}: respuesta vacía")
            continue
        # rows[0] = nombres, rows[1] = unidades, rows[2] = guiones
        data = list(csv.reader(io.StringIO("\n".join(rows[3:])), delimiter="\t"))

        def col(i):
            return np.array([float(d[i]) if len(d) > i and d[i].strip() else np.nan for d in data])

        tbl = Table({"ra": col(0), "dec": col(1), "V": col(2), "eV": col(3), "B": col(4), "eB": col(5)})
        tbl = tbl[np.isfinite(tbl["B"]) & np.isfinite(tbl["V"])]
        tbl.meta["servidor"] = server
        log.info("APASS: %d estrellas con B y V en %.0f′ (%s)", len(tbl), radius_arcmin, server)
        return tbl
    raise ConnectionError("No se pudo consultar APASS en VizieR:\n  " + "\n  ".join(errors))


def match_field_stars(cmd, ref, target, exclude_arcmin=6.0, vmin=11.0, vmax=16.0, max_sep_arcsec=2.0,
                      max_err=0.1):
    """Cruza el catálogo B−V calibrado con el de referencia y elige las estrellas útiles.

    Se quedan las estrellas no saturadas, con errores < ``max_err`` en ambos catálogos,
    ``vmin`` < V < ``vmax`` (APASS) y a más de ``exclude_arcmin`` del objeto (la galaxia contamina
    la fotometría). Devuelve una tabla con las magnitudes de ambos catálogos y las diferencias.
    """
    pos = SkyCoord(np.asarray(cmd["ra"]) * u.deg, np.asarray(cmd["dec"]) * u.deg)
    cat = SkyCoord(np.asarray(ref["ra"]) * u.deg, np.asarray(ref["dec"]) * u.deg)
    idx, sep, _ = pos.match_to_catalog_sky(cat)
    r = ref[idx]
    ok = ((sep.arcsec < max_sep_arcsec) & ~np.asarray(cmd["saturated"], bool)
          & np.isfinite(cmd["B"]) & np.isfinite(cmd["V"])
          & (np.asarray(cmd["B_err"]) < max_err) & (np.asarray(cmd["V_err"]) < max_err)
          & (np.asarray(r["eB"]) < max_err) & (np.asarray(r["eV"]) < max_err)
          & (np.asarray(r["V"]) > vmin) & (np.asarray(r["V"]) < vmax)
          & (pos.separation(target).arcmin > exclude_arcmin))
    out = Table({"x": cmd["x"][ok], "y": cmd["y"][ok], "ra": cmd["ra"][ok], "dec": cmd["dec"][ok],
                 "sep_arcsec": sep.arcsec[ok],
                 "B_cat": r["B"][ok], "V_cat": r["V"][ok], "B_V_cat": r["B"][ok] - r["V"][ok],
                 "B": cmd["B"][ok], "V": cmd["V"][ok], "B_V": cmd["B_V"][ok]})
    out["dB"] = out["B_cat"] - out["B"]
    out["dV"] = out["V_cat"] - out["V"]
    return out


def offsets(stars, sigma=2.5):
    """Δ por filtro: mediana con sigma-clipping, rms y error de la media (rms/√N)."""
    res = {}
    for filt in ("B", "V"):
        d = np.asarray(stars[f"d{filt}"], float)
        if len(d) < 3:
            res[filt] = {"delta": np.nan, "rms": np.nan, "err": np.nan, "n": len(d)}
            continue
        _, med, std = sigma_clipped_stats(d, sigma=sigma, maxiters=5)
        n = int(np.sum(np.abs(d - med) < sigma * max(std, 1e-3)))
        res[filt] = {"delta": float(med), "rms": float(std), "err": float(std / np.sqrt(n)), "n": n}
    return res


def corrected_zeropoints(zps, delta):
    """Zeropoints de la noche corregidos por transparencia (ver la fórmula del módulo)."""
    out = {f: dict(v) for f, v in zps.items()}
    dv = delta["V"]["delta"]
    if "V" in out and np.isfinite(dv):
        out["V"]["zp"] += dv
        out["V"]["zp_err"] = float(np.hypot(out["V"]["zp_err"], delta["V"]["err"]))
    if "B" in out and np.isfinite(delta["B"]["delta"]):
        c = out["B"].get("c", 0.0)
        db = delta["B"]["delta"]
        shift = db * (1 - c) + c * (dv if np.isfinite(dv) else db)
        out["B"]["zp"] += shift
        out["B"]["zp_err"] = float(np.hypot(out["B"]["zp_err"], delta["B"]["err"]))
    return out
