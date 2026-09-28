"""Detección de fuentes y generación de catálogos con photutils.

La tarea pide (al menos) dos configuraciones que difieran en el *deblending*:

- ``agresiva``: fondo local fino (caja chica) más un filtro de mediana que
  quita las estructuras mayores a unos pocos FWHM (la nebulosa), filtro del
  tamaño del seeing, umbral moderado, muchos niveles y contraste muy bajo.
  Separa estrellas vecinas dentro del cúmulo. Las
  magnitudes se miden con fotometría de apertura en las posiciones detectadas.
- ``extendida``: fondo con caja grande (la nebulosa **no** se resta),
  filtro ancho, umbral bajo, área mínima grande y contraste alto. Detecta la
  emisión difusa como pocas regiones grandes, sin fragmentarla en "fuentes".

Los parámetros por defecto están en ``DEFAULT_CONFIGS`` y se pueden cambiar
desde el YAML (sección ``deteccion``).
"""

import warnings

import numpy as np
from astropy.stats import mad_std
from astropy.table import Table
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import SourceCatalog
from photutils.utils import calc_total_error
from scipy import ndimage
from scipy.spatial import cKDTree

from . import compat

FWHM_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))

DEFAULT_CONFIGS = {
    "agresiva": {
        "box_size": 64,            # px, caja del fondo local
        "escala_nebulosa_fwhm": 4, # resta estructuras > 4 FWHM (la nebulosa), con las fuentes enmascaradas
        "n_sigma": 3.0,            # umbral en sigmas de la imagen filtrada
        "n_pixels": 10,            # área mínima conectada (px)
        "kernel_fwhm_factor": 1.0, # FWHM del filtro = factor * seeing
        "deblend": True,
        "n_levels": 64,
        "contrast": 0.0005,
        "flux": "aperture",        # apertura chica + fondo local + corrección de apertura
        "r_apertura_fwhm": 1.0,    # radio de la apertura en unidades de FWHM
    },
    "extendida": {
        "box_size": 512,
        "n_sigma": 1.5,
        "n_pixels": 400,
        "kernel_fwhm_factor": 3.0,
        "deblend": True,
        "n_levels": 16,
        "contrast": 0.2,
        "flux": "segment",
    },
}


def subtract_background(data, box_size=64, filter_size=5):
    """Estima el fondo con ``Background2D``. Devuelve ``(data - fondo, fondo, rms)``."""
    ny, nx = data.shape
    box = int(min(box_size, ny // 2, nx // 2))
    coverage = ~np.isfinite(data)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bkg = Background2D(data, box, coverage_mask=coverage, filter_size=filter_size,
                           bkg_estimator=MedianBackground(), exclude_percentile=50.0, fill_value=np.nan)
    background = np.asarray(bkg.background, dtype=np.float32)
    rms = np.asarray(bkg.background_rms, dtype=np.float32)
    return (data - background).astype(np.float32), background, rms


def large_scale(data, scale_px):
    """Estructura de escala mayor que ``scale_px`` (nebulosa), por filtro de mediana.

    Para que sea rápido en 4096x4096, la mediana se calcula sobre la imagen
    binneada y luego se interpola de vuelta al tamaño original.
    """
    ny, nx = data.shape
    b = max(1, int(scale_px // 6))
    my, mx = ny // b * b, nx // b * b
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        binned = np.nanmedian(data[:my, :mx].reshape(my // b, b, mx // b, b), axis=(1, 3))
    bad = ~np.isfinite(binned)
    if bad.all():
        return np.zeros_like(data, dtype=np.float32)
    if bad.any():
        # Celdas sin datos (bordes o fuentes enmascaradas): valor de la celda válida más cercana
        idx = ndimage.distance_transform_edt(bad, return_distances=False, return_indices=True)
        binned = binned[tuple(idx)]
    size = max(3, int(round(scale_px / b)) | 1)
    smooth = ndimage.median_filter(binned, size=size, mode="nearest")
    smooth = ndimage.gaussian_filter(smooth, size / 4)
    zoomed = ndimage.zoom(smooth, (ny / smooth.shape[0], nx / smooth.shape[1]), order=1)
    return zoomed[:ny, :nx].astype(np.float32)


def detect(data_sub, rms, fwhm_px, params):
    """Filtra con una gaussiana, umbraliza y (opcionalmente) separa fuentes.

    El umbral se escala por el ruido de la imagen **filtrada** (que es menor
    que el de la original), así ``n_sigma`` es la significancia real de la
    detección. Devuelve ``(segmentación o None, imagen filtrada)``.
    """
    good = np.isfinite(data_sub)
    sigma_k = params.get("kernel_fwhm_factor", 1.0) * fwhm_px / FWHM_SIGMA
    filled = np.where(good, data_sub, 0.0)
    conv = ndimage.gaussian_filter(filled, sigma_k).astype(np.float32)
    ratio = mad_std(conv[good]) / mad_std(filled[good])
    threshold = params["n_sigma"] * np.where(np.isfinite(rms), rms, np.nanmedian(rms)) * ratio
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        segm = compat.detect_sources(conv, threshold, params["n_pixels"], mask=~good)
        if segm is not None and params.get("deblend", False):
            segm = compat.deblend_sources(conv, segm, params["n_pixels"],
                                          n_levels=params.get("n_levels", 32),
                                          contrast=params.get("contrast", 0.001))
    return segm, conv


def source_properties(data_sub, segm, conv, error=None, kron=False):
    """Tabla con las propiedades básicas de cada segmento."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cat = SourceCatalog(data_sub, segm, convolved_data=conv, error=error, mask=~np.isfinite(data_sub))
        tbl = Table()
        tbl["label"] = np.asarray(cat.label)
        tbl["x"] = compat.cat_value(cat, "x_centroid", "xcentroid")
        tbl["y"] = compat.cat_value(cat, "y_centroid", "ycentroid")
        tbl["area"] = compat.cat_value(cat, "area")
        tbl["flux"] = compat.cat_value(cat, "segment_flux")
        tbl["flux_err"] = compat.cat_value(cat, "segment_flux_err") if error is not None else np.nan
        tbl["peak"] = compat.cat_value(cat, "max_value")
        tbl["elongation"] = compat.cat_value(cat, "elongation")
        if kron:
            tbl["kron_flux"] = compat.cat_value(cat, "kron_flux")
            tbl["kron_flux_err"] = compat.cat_value(cat, "kron_flux_err") if error is not None else np.nan
    return tbl


def build_catalog(data, fwhm_px, params, gain_eff=None, sat_mask=None, wcs=None):
    """Detecta fuentes y arma el catálogo de una imagen apilada (ADU/s).

    ``gain_eff`` = ganancia (e-/ADU) x tiempo total de exposición; se usa para
    el error de Poisson de las fuentes.

    Devuelve ``(tabla, segmentación, data_sub, fondo, rms)``. La columna
    ``flux`` es el flujo Kron (``params["flux"] == "kron"``) o el flujo total
    del segmento (``"segment"``, para estructuras extendidas). Con
    ``"aperture"`` queda el flujo del segmento y el pipeline lo reemplaza por
    fotometría de apertura (``photometry.aperture_photometry_list``), que es
    la más confiable para estrellas en un campo con nebulosa.
    """
    data_sub, background, rms = subtract_background(data, box_size=params.get("box_size", 64))
    scale = params.get("escala_nebulosa_fwhm")
    if scale:
        # 1) Detección preliminar sobre la imagen sin estructuras mayores a `scale` FWHM
        segm0, _ = detect(data_sub - large_scale(data_sub, scale * fwhm_px), rms, fwhm_px,
                          {**params, "deblend": False})
        # 2) Mapa de la nebulosa con las fuentes enmascaradas, para que las estrellas no
        #    "inflen" la mediana y se resten a sí mismas parte de su flujo
        masked = data_sub.copy()
        if segm0 is not None:
            masked[ndimage.binary_dilation(segm0.data > 0, iterations=int(np.ceil(fwhm_px)))] = np.nan
        nebula = large_scale(masked, scale * fwhm_px)
        data_sub = data_sub - nebula
        background = background + nebula
    segm, conv = detect(data_sub, rms, fwhm_px, params)
    if segm is None:
        return Table(), None, data_sub, background, rms
    error = None
    if gain_eff:
        error = calc_total_error(np.where(np.isfinite(data_sub), data_sub, 0.0), rms, gain_eff)
    use_kron = params.get("flux", "kron") == "kron"
    tbl = source_properties(data_sub, segm, conv, error=error, kron=use_kron)
    if use_kron:
        kron_ok = np.isfinite(tbl["kron_flux"]) & (tbl["kron_flux"] > 0)
        tbl["flux"] = np.where(kron_ok, tbl["kron_flux"], tbl["flux"])
        tbl["flux_err"] = np.where(kron_ok, tbl["kron_flux_err"], tbl["flux_err"])
        tbl.remove_columns(["kron_flux", "kron_flux_err"])

    tbl["saturated"] = False
    if sat_mask is not None and sat_mask.any():
        bad = np.unique(segm.data[sat_mask & (segm.data > 0)])
        tbl["saturated"] = np.isin(tbl["label"], bad)
    if wcs is not None:
        ra, dec = wcs.pixel_to_world_values(tbl["x"], tbl["y"])
        tbl["ra"], tbl["dec"] = ra, dec
    return tbl, segm, data_sub, background, rms


def calibrate(tbl, zp, zp_err, k, airmass):
    """Agrega magnitudes calibradas: m = ZP - k X - 2.5 log10(F [ADU/s])."""
    with np.errstate(invalid="ignore", divide="ignore"):
        flux = np.asarray(tbl["flux"], dtype=float)
        tbl["mag"] = zp - k * airmass - 2.5 * np.log10(flux)
        ferr = np.asarray(tbl["flux_err"], dtype=float)
        tbl["mag_err"] = np.hypot(2.5 / np.log(10) * ferr / flux, zp_err)
        bad = ~(flux > 0)
        tbl["mag"][bad] = np.nan
        tbl["mag_err"][bad] = np.nan
    return tbl


def crossmatch(tbl_b, tbl_v, tol_px):
    """Une los catálogos B y V (en la misma grilla de píxeles) por posición.

    Cada fuente de B se asocia a la fuente de V más cercana dentro de
    ``tol_px``; si dos fuentes de B caen en la misma de V se queda la más
    cercana. Devuelve una tabla con B, V y B-V.
    """
    pts_v = np.column_stack([tbl_v["x"], tbl_v["y"]])
    pts_b = np.column_stack([tbl_b["x"], tbl_b["y"]])
    dist, idx = cKDTree(pts_v).query(pts_b, distance_upper_bound=tol_px)
    matched = np.isfinite(dist)
    order = np.argsort(dist)
    used, rows = set(), []
    for i in order:
        if not matched[i] or idx[i] in used:
            continue
        used.add(idx[i])
        rows.append((i, idx[i], dist[i]))
    ib = np.array([r[0] for r in rows], dtype=int)
    iv = np.array([r[1] for r in rows], dtype=int)
    out = Table()
    out["x"] = np.asarray(tbl_v["x"])[iv]
    out["y"] = np.asarray(tbl_v["y"])[iv]
    for col in ("ra", "dec"):
        if col in tbl_v.colnames:
            out[col] = np.asarray(tbl_v[col])[iv]
    out["sep_px"] = [r[2] for r in rows]
    out["B"] = np.asarray(tbl_b["mag"])[ib]
    out["B_err"] = np.asarray(tbl_b["mag_err"])[ib]
    out["V"] = np.asarray(tbl_v["mag"])[iv]
    out["V_err"] = np.asarray(tbl_v["mag_err"])[iv]
    out["B_V"] = out["B"] - out["V"]
    out["B_V_err"] = np.hypot(out["B_err"], out["V_err"])
    sat_b = np.asarray(tbl_b["saturated"])[ib] if "saturated" in tbl_b.colnames else False
    sat_v = np.asarray(tbl_v["saturated"])[iv] if "saturated" in tbl_v.colnames else False
    out["saturated"] = np.logical_or(sat_b, sat_v)
    return out


def limiting_magnitude(rms, fwhm_px, zp_eff, n_sigma=5.0):
    """Magnitud límite para una fuente puntual detectada a ``n_sigma``.

    Apertura de radio = FWHM (contiene el 93.75 % del flujo de una gaussiana).
    ``zp_eff`` = ZP - k X. ``rms`` debe venir corregido por NOISECOR (el
    alineamiento interpola y correlaciona el ruido del apilado).
    """
    area = np.pi * fwhm_px ** 2
    flux = n_sigma * rms * np.sqrt(area) / 0.9375
    return float(zp_eff - 2.5 * np.log10(flux))
