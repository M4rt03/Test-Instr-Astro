"""Estimación del seeing a partir de estrellas de la imagen de ciencia.

A diferencia de medir una sola estrella estándar (que además puede estar
saturada), aquí se ajusta una gaussiana 2D elíptica a decenas de estrellas
aisladas y no saturadas de cada exposición. Eso da:

- el seeing típico (mediana del FWHM) y su dispersión,
- el seeing de cada exposición (para comparar noches y filtros),
- la variación espacial de la PSF en el campo (aberraciones, foco, tracking).
"""

import warnings

import numpy as np
from astropy.stats import mad_std, sigma_clip
from astropy.table import Table
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

from . import detection

FWHM_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))  # 2.3548

# Longitud de onda efectiva de los filtros Johnson (nm)
LAMBDA_NM = {"U": 365.0, "B": 445.0, "V": 551.0, "R": 658.0, "I": 806.0}


def _gauss2d(p, y, x):
    amp, x0, y0, sx, sy, theta, offset = p
    ct, st = np.cos(theta), np.sin(theta)
    xr = (x - x0) * ct + (y - y0) * st
    yr = -(x - x0) * st + (y - y0) * ct
    return offset + amp * np.exp(-0.5 * ((xr / sx) ** 2 + (yr / sy) ** 2))


def fit_gaussian(cutout, sigma_guess):
    """Ajusta una gaussiana 2D elíptica rotada + fondo constante a un recorte.

    Devuelve un dict con centro (en coordenadas del recorte), FWHM en cada eje,
    FWHM medio geométrico y elipticidad, o ``None`` si el ajuste falla.
    """
    ny, nx = cutout.shape
    y, x = np.mgrid[:ny, :nx]
    good = np.isfinite(cutout)
    if good.sum() < 0.8 * cutout.size:
        return None
    data = cutout[good]
    yy, xx = y[good], x[good]
    offset0 = np.median(np.concatenate([cutout[0][np.isfinite(cutout[0])], cutout[-1][np.isfinite(cutout[-1])]]))
    amp0 = np.nanmax(cutout) - offset0
    p0 = [amp0, (nx - 1) / 2, (ny - 1) / 2, sigma_guess, sigma_guess, 0.0, offset0]
    lower = [0, 0, 0, 0.3, 0.3, -np.pi, -np.inf]
    upper = [np.inf, nx - 1, ny - 1, nx, ny, np.pi, np.inf]
    try:
        res = least_squares(lambda p: _gauss2d(p, yy, xx) - data, p0, bounds=(lower, upper), max_nfev=200)
    except ValueError:
        return None
    if not res.success:
        return None
    amp, x0, y0, sx, sy, theta, offset = res.x
    fx, fy = FWHM_SIGMA * sx, FWHM_SIGMA * sy
    major, minor = max(fx, fy), min(fx, fy)
    return {
        "x0": x0, "y0": y0, "amp": amp, "offset": offset,
        "fwhm_major": major, "fwhm_minor": minor,
        "fwhm": np.sqrt(major * minor),
        "ellipticity": 1 - minor / major,
    }


def measure_psf(data, sat_mask=None, fwhm_guess_px=8.0, n_sigma=10.0, max_stars=150,
                isolation_factor=4.0, box_size=64):
    """Mide el FWHM de las estrellas aisladas y no saturadas de una imagen.

    Devuelve una tabla con una fila por estrella (``x``, ``y``, ``fwhm_px``,
    ``fwhm_major``, ``fwhm_minor``, ``ellipticity``, ``peak``).
    """
    data_sub, _, rms = detection.subtract_background(data, box_size=box_size)
    data_sub = data_sub - detection.large_scale(data_sub, 5 * fwhm_guess_px)  # quita la nebulosa
    params = {"n_sigma": n_sigma, "n_pixels": 8, "kernel_fwhm_factor": 1.0, "deblend": False}
    segm, conv = detection.detect(data_sub, rms, fwhm_guess_px, params)
    if segm is None:
        return Table(names=("x", "y", "fwhm_px", "fwhm_major", "fwhm_minor", "ellipticity", "peak"))
    cat = detection.source_properties(data_sub, segm, conv)

    x, y, flux = cat["x"], cat["y"], cat["flux"]
    ok = np.isfinite(x) & np.isfinite(y) & (flux > 0)
    # Excluye estrellas que tocan píxeles saturados
    if sat_mask is not None and sat_mask.any():
        bad_labels = np.unique(segm.data[sat_mask & (segm.data > 0)])
        ok &= ~np.isin(cat["label"], bad_labels)
    # Aisladas: sin otra fuente detectada a menos de isolation_factor * FWHM
    pts = np.column_stack([x, y])
    if len(pts) > 1:
        dist, _ = cKDTree(pts).query(pts, k=2)
        ok &= dist[:, 1] > isolation_factor * fwhm_guess_px

    half = int(max(6, round(2.5 * fwhm_guess_px)))
    ny, nx = data.shape
    ok &= (x > half) & (x < nx - half - 1) & (y > half) & (y < ny - half - 1)

    idx = np.where(ok)[0]
    idx = idx[np.argsort(flux[idx])[::-1]][:max_stars]

    rows = []
    for i in idx:
        xi, yi = int(round(x[i])), int(round(y[i]))
        cut = data_sub[yi - half:yi + half + 1, xi - half:xi + half + 1]
        fit = fit_gaussian(cut, fwhm_guess_px / FWHM_SIGMA)
        if fit is None:
            continue
        rows.append((xi - half + fit["x0"], yi - half + fit["y0"], fit["fwhm"], fit["fwhm_major"],
                     fit["fwhm_minor"], fit["ellipticity"], fit["amp"]))
    tbl = Table(rows=rows, names=("x", "y", "fwhm_px", "fwhm_major", "fwhm_minor", "ellipticity", "peak"))
    if len(tbl) == 0:
        return tbl
    # Descarta ajustes absurdos (rayos cósmicos, galaxias, blends) y outliers
    sane = (tbl["fwhm_px"] > 1.0) & (tbl["fwhm_px"] < 5 * fwhm_guess_px) & (tbl["ellipticity"] < 0.5)
    tbl = tbl[sane]
    if len(tbl) > 5:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            clipped = sigma_clip(np.asarray(tbl["fwhm_px"]), sigma=3, maxiters=5)
        tbl = tbl[~clipped.mask]
    return tbl


def summarize(tbl, pixscale):
    """Resumen robusto del FWHM de una tabla de ``measure_psf``."""
    if len(tbl) == 0:
        return {"n_estrellas": 0, "fwhm_px": np.nan, "seeing_arcsec": np.nan,
                "dispersion_arcsec": np.nan, "elipticidad": np.nan}
    fwhm = np.asarray(tbl["fwhm_px"])
    return {
        "n_estrellas": len(tbl),
        "fwhm_px": float(np.median(fwhm)),
        "seeing_arcsec": float(np.median(fwhm) * pixscale),
        "dispersion_arcsec": float(mad_std(fwhm) * pixscale),
        "elipticidad": float(np.median(tbl["ellipticity"])),
    }


def diffraction_limit(filt, diameter_m):
    """Límite de difracción en arcsec: radio del primer anillo oscuro (1.22 λ/D) y FWHM de Airy (1.03 λ/D)."""
    lam = LAMBDA_NM.get(filt, 550.0) * 1e-9
    rad_to_arcsec = 180 / np.pi * 3600
    return {
        "lambda_nm": lam * 1e9,
        "rayleigh_arcsec": 1.22 * lam / diameter_m * rad_to_arcsec,
        "fwhm_airy_arcsec": 1.029 * lam / diameter_m * rad_to_arcsec,
    }


def spatial_grid(tbl, shape, n=3, min_stars=3):
    """Mediana del FWHM en una grilla de ``n`` x ``n`` celdas sobre el detector.

    Si la PSF varía de forma importante (por ejemplo, más ancha en las
    esquinas por curvatura de campo o coma), se ve aquí.
    """
    ny, nx = shape
    grid = np.full((n, n), np.nan)
    counts = np.zeros((n, n), dtype=int)
    if len(tbl) == 0:
        return grid, counts
    ix = np.clip((np.asarray(tbl["x"]) / nx * n).astype(int), 0, n - 1)
    iy = np.clip((np.asarray(tbl["y"]) / ny * n).astype(int), 0, n - 1)
    for j in range(n):
        for i in range(n):
            sel = (ix == i) & (iy == j)
            counts[j, i] = sel.sum()
            if counts[j, i] >= min_stars:
                grid[j, i] = np.median(np.asarray(tbl["fwhm_px"])[sel])
    return grid, counts


def radial_trend(tbl, shape):
    """Pendiente del FWHM (px) contra la distancia al centro (px por cada 1000 px).

    Una pendiente positiva y significativa indica que la PSF empeora hacia
    los bordes del campo.
    """
    if len(tbl) < 5:
        return np.nan, np.nan
    ny, nx = shape
    r = np.hypot(np.asarray(tbl["x"]) - nx / 2, np.asarray(tbl["y"]) - ny / 2) / 1000.0
    coef, cov = np.polyfit(r, np.asarray(tbl["fwhm_px"]), 1, cov=True)
    return float(coef[0]), float(np.sqrt(cov[0, 0]))
