"""Calibración fotométrica con estrellas estándar.

Relación usada (sin término de color, porque con dos estrellas no se puede ajustar):

    m_std = -2.5 log10(F) + ZP - k X

- F  : flujo de la estrella en ADU/s (apertura grande, fondo local restado).
- ZP : zeropoint del filtro (lo que se mide).
- k  : coeficiente de extinción (mag/masa de aire). Con dos estrellas a
       masas de aire parecidas no se puede ajustar bien, así que se usa un
       valor típico de Cerro Tololo (configurable en el YAML).
- X  : masa de aire (del header o calculada a partir de la hora y posición).

Importante: las magnitudes del catálogo Tycho (BT, VT) **no** son Johnson
B, V; hay que convertirlas (``tycho_to_johnson``).
"""

import warnings

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.table import Table
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry
from scipy import ndimage

from .io import obs_time
from .seeing import FWHM_SIGMA, fit_gaussian

LN10 = np.log(10.0)


def tycho_to_johnson(bt, vt):
    """Convierte Tycho (BT, VT) a Johnson (B, V).

    Transformación de ESA (1997), The Hipparcos and Tycho Catalogues, §1.3:
    V = VT - 0.090 (BT - VT),  B - V = 0.850 (BT - VT);
    válida para -0.2 < BT - VT < 1.8.
    """
    bt_vt = bt - vt
    if not -0.2 < bt_vt < 1.8:
        warnings.warn(f"BT-VT = {bt_vt:.2f} fuera del rango de validez de la conversión")
    v = vt - 0.090 * bt_vt
    return v + 0.850 * bt_vt, v


def standard_magnitudes(std_cfg):
    """Magnitudes Johnson de una estándar desde el YAML (claves B/V o BT/VT)."""
    if "B" in std_cfg.get("magnitudes", {}):
        mags = std_cfg["magnitudes"]
        return {"B": float(mags["B"]), "V": float(mags["V"])}
    mags = std_cfg["magnitudes"]
    b, v = tycho_to_johnson(float(mags["BT"]), float(mags["VT"]))
    return {"B": b, "V": v}


def site_location(site_cfg):
    return EarthLocation(lat=site_cfg["latitud"] * u.deg, lon=site_cfg["longitud"] * u.deg,
                         height=site_cfg["altitud"] * u.m)


def airmass(header, location=None, coord=None):
    """Masa de aire: la del header si existe; si no, se calcula (sec z)."""
    value = header.get("AIRMASS")
    if value not in (None, "") and np.isfinite(float(value)) and float(value) >= 1:
        return float(value)
    if location is None or coord is None:
        raise ValueError("El header no trae AIRMASS; hay que entregar sitio y coordenadas")
    altaz = coord.transform_to(AltAz(obstime=obs_time(header), location=location))
    return float(altaz.secz)


def aperture_flux(data, x, y, r_ap, r_in, r_out):
    """Suma en una apertura circular menos el fondo local (mediana clipeada del anillo)."""
    aper = CircularAperture((x, y), r=r_ap)
    annulus = CircularAnnulus((x, y), r_in=r_in, r_out=r_out)
    mask = ~np.isfinite(data)
    phot = aperture_photometry(data, aper, mask=mask)
    ann_vals = annulus.to_mask(method="center").get_values(data)
    ann_vals = ann_vals[np.isfinite(ann_vals)]
    _, sky, sky_std = sigma_clipped_stats(ann_vals, sigma=3.0)
    flux = float(phot["aperture_sum"][0]) - sky * aper.area
    return flux, float(sky), float(sky_std), float(aper.area)


def aperture_photometry_list(data, x, y, r_ap, r_in, r_out, gain_eff=None, noise_corr=1.0):
    """Fotometría de apertura para muchas fuentes a la vez, con fondo local de cada una.

    El fondo es la mediana clipeada (3σ) de un anillo alrededor de cada
    fuente, así se adapta a la nebulosa. Devuelve ``(flujo, error, fondo)``.
    ``noise_corr`` (header NOISECOR del apilado) corrige el ruido por píxel,
    que la interpolación del alineamiento hace parecer menor de lo que es.
    El error incluye el ruido de Poisson de la fuente (si se da ``gain_eff``),
    el ruido del cielo en la apertura, la incertidumbre estadística del fondo
    y un término por estructura del fondo: la diferencia entre la mitad
    interna y la externa del anillo. Sobre una nebulosa con curvatura esa
    diferencia es grande y el error crece, como corresponde.
    """
    from astropy.stats import SigmaClip
    from photutils.aperture import ApertureStats

    pos = np.column_stack([x, y])
    aper = CircularAperture(pos, r=r_ap)
    annulus = CircularAnnulus(pos, r_in=r_in, r_out=r_out)
    mask = ~np.isfinite(data)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        phot = aperture_photometry(data, aper, mask=mask)
        ann = ApertureStats(data, annulus, mask=mask, sigma_clip=SigmaClip(sigma=3.0, maxiters=5))
        r_mid = np.sqrt(0.5 * (r_in ** 2 + r_out ** 2))  # divide el anillo en dos de igual área
        clip = SigmaClip(sigma=3.0, maxiters=5)
        inner = ApertureStats(data, CircularAnnulus(pos, r_in=r_in, r_out=r_mid), mask=mask, sigma_clip=clip)
        outer = ApertureStats(data, CircularAnnulus(pos, r_in=r_mid, r_out=r_out), mask=mask, sigma_clip=clip)
    sky = np.asarray(ann.median, dtype=float)
    sky_std = np.asarray(ann.std, dtype=float) * noise_corr
    n_sky = np.asarray(getattr(ann.center_aper_area, "value", ann.center_aper_area), dtype=float)
    area = aper.area
    flux = np.asarray(phot["aperture_sum"], dtype=float) - sky * area
    stat_err = sky_std * np.sqrt(2.0 / np.maximum(n_sky, 1))  # ruido esperado de la diferencia
    structure = np.asarray(inner.median, dtype=float) - np.asarray(outer.median, dtype=float)
    struct_var = np.maximum(structure ** 2 - stat_err ** 2, 0.0)  # solo el exceso sobre el ruido
    var = area * sky_std ** 2 + area ** 2 * sky_std ** 2 / np.maximum(n_sky, 1) + area ** 2 * struct_var
    if gain_eff:
        var = var + np.maximum(flux, 0) / gain_eff
    return flux, np.sqrt(var), sky


def _growth_ratio(stamp, c, r_small, r_big):
    small = aperture_photometry(stamp, CircularAperture((c, c), r=r_small))["aperture_sum"][0]
    big = aperture_photometry(stamp, CircularAperture((c, c), r=r_big))["aperture_sum"][0]
    return float(big / small)


def _combine_stamps(stamps, weights):
    """Promedio pesado píxel a píxel, rechazando valores a más de 3σ (vecinas, rayos cósmicos)."""
    clipped = sigma_clip(stamps, sigma=3.0, maxiters=3, axis=0, cenfunc="median", stdfunc="mad_std")
    w = np.where(clipped.mask, 0.0, weights[:, None, None])
    return np.sum(np.where(clipped.mask, 0.0, stamps) * w, axis=0) / np.maximum(np.sum(w, axis=0), 1e-30)


def aperture_correction(data, x, y, r_small, r_big, r_in, r_out, max_stars=40, min_snr=20.0):
    """Factor F(r_big) / F(r_small) a partir de una PSF empírica (curva de crecimiento).

    El catálogo se mide con apertura chica (menos contaminación de vecinas y
    menos ruido) y se lleva al flujo total, el mismo con que se calibró el ZP
    en las estándares, multiplicando por este factor.

    Para cada estrella brillante (S/N > ``min_snr`` en la apertura chica) se
    toma un recorte centrado con precisión sub-píxel, se le resta el fondo del
    anillo y se normaliza por su flujo en ``r_small``. La PSF empírica es el
    promedio de los recortes pesado por (S/N)², con rechazo 3σ píxel a píxel:
    una vecina aparece solo en uno o dos recortes y se rechaza, así no hace
    falta que las estrellas estén perfectamente aisladas. El error se estima
    con jackknife (dejando fuera una estrella a la vez).

    Devuelve ``(factor, error del factor, n estrellas usadas)``.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if len(x) < 3:
        return 1.0, np.nan, len(x)
    f_small, e_small, _ = aperture_photometry_list(data, x, y, r_small, r_in, r_out)
    snr = np.nan_to_num(f_small / e_small, nan=-np.inf)
    order = [i for i in np.argsort(snr)[::-1][:max_stars] if snr[i] > min_snr]

    half = int(np.ceil(r_out)) + 2
    c = float(half)
    yy, xx = np.mgrid[:2 * half + 1, :2 * half + 1]
    rr = np.hypot(xx - c, yy - c)
    ring = (rr >= r_in) & (rr <= r_out)
    stamps, weights = [], []
    for i in order:
        xi, yi = int(round(x[i])), int(round(y[i]))
        if xi - half < 0 or yi - half < 0 or xi + half + 1 > data.shape[1] or yi + half + 1 > data.shape[0]:
            continue
        cut = np.array(data[yi - half:yi + half + 1, xi - half:xi + half + 1], dtype=float)
        if not np.all(np.isfinite(cut)):
            continue
        cut = ndimage.shift(cut, (yi - y[i], xi - x[i]), order=3, mode="nearest")
        _, sky, _ = sigma_clipped_stats(cut[ring], sigma=3.0)
        cut -= sky
        norm = aperture_photometry(cut, CircularAperture((c, c), r=r_small))["aperture_sum"][0]
        if norm > 0:
            stamps.append(cut / norm)
            weights.append(snr[i] ** 2)
    if len(stamps) < 3:
        return 1.0, np.nan, len(stamps)
    stamps, weights = np.array(stamps), np.array(weights)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        factor = _growth_ratio(_combine_stamps(stamps, weights), c, r_small, r_big)
        n = len(stamps)
        keep = ~np.eye(n, dtype=bool)
        jack = np.array([_growth_ratio(_combine_stamps(stamps[k], weights[k]), c, r_small, r_big) for k in keep])
    return factor, float(np.sqrt((n - 1) / n * np.sum((jack - jack.mean()) ** 2))), n


def measure_standard(data, sat_mask=None, fwhm_guess_px=8.0, aperture_factor=3.0):
    """Localiza la estrella más brillante (la estándar) y mide su flujo.

    - Posición: máximo de la imagen suavizada (evita rayos cósmicos) y
      refinamiento con un ajuste gaussiano, que además da el FWHM.
    - Flujo: apertura de radio ``aperture_factor`` x FWHM (>99 % del flujo de
      una PSF gaussiana) y anillo de fondo entre ``aperture_factor`` + 0.5 y
      ``aperture_factor`` + 2 FWHM (el mismo que usa el catálogo).
    """
    good = np.isfinite(data)
    filled = np.where(good, data, np.nanmedian(data))
    smooth = ndimage.gaussian_filter(filled, fwhm_guess_px / FWHM_SIGMA)
    yi, xi = np.unravel_index(np.argmax(smooth), smooth.shape)

    half = int(max(8, round(3 * fwhm_guess_px)))
    y0, x0 = max(yi - half, 0), max(xi - half, 0)
    cut = data[y0:yi + half + 1, x0:xi + half + 1]
    fit = fit_gaussian(cut, fwhm_guess_px / FWHM_SIGMA)
    if fit is not None:
        x, y, fwhm = x0 + fit["x0"], y0 + fit["y0"], fit["fwhm"]
    else:
        x, y, fwhm = float(xi), float(yi), fwhm_guess_px
    r_ap = aperture_factor * fwhm
    flux, sky, sky_std, area = aperture_flux(data, x, y, r_ap, (aperture_factor + 0.5) * fwhm,
                                             (aperture_factor + 2.0) * fwhm)

    saturated = False
    if sat_mask is not None:
        yy, xx = np.ogrid[:data.shape[0], :data.shape[1]]
        near = (xx - x) ** 2 + (yy - y) ** 2 <= (2 * fwhm) ** 2
        saturated = bool(np.any(sat_mask & near))
    return {"x": x, "y": y, "fwhm_px": fwhm, "r_ap": r_ap, "flux": flux,
            "sky": sky, "sky_std": sky_std, "saturated": saturated}


def zeropoint(flux, m_std, k, x_airmass):
    """ZP = m_std + 2.5 log10(F) + k X."""
    return m_std + 2.5 * np.log10(flux) + k * x_airmass


def summarize_zeropoints(tbl):
    """ZP por filtro: media clipeada, desviación estándar y error de la media (σ/√N).

    Se excluyen las mediciones de estrellas saturadas.
    """
    rows = []
    for filt in sorted(set(tbl["filtro"])):
        sel = (tbl["filtro"] == filt) & ~np.asarray(tbl["saturada"], dtype=bool) & np.isfinite(tbl["zp"])
        zps = np.asarray(tbl["zp"][sel])
        if len(zps) == 0:
            rows.append((filt, np.nan, np.nan, np.nan, 0))
            continue
        mean, _, std = sigma_clipped_stats(zps, sigma=3.0)
        n = len(zps)
        rows.append((filt, float(mean), float(std), float(std / np.sqrt(n)) if n > 1 else np.nan, n))
    return Table(rows=rows, names=("filtro", "zp", "zp_std", "zp_err", "n"))
