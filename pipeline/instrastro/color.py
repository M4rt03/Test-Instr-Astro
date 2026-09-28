"""Imagen en color y mapa de color B-V.

Con solo dos filtros se asigna R = V, B = B y un canal verde sintético
G = (B + V) / 2, como sugiere el enunciado. Antes de combinar se resta el
fondo de cielo a cada canal y se los escala al mismo nivel, para que el color
refleje diferencias reales entre filtros y no diferencias de fondo o de
tiempo de exposición.
"""

import numpy as np
from astropy.stats import sigma_clipped_stats
from astropy.visualization import make_lupton_rgb
from scipy import ndimage


def _sky_subtract(img):
    good = np.isfinite(img)
    _, median, _ = sigma_clipped_stats(img[good], sigma=3.0, maxiters=5)
    return np.where(good, img - median, 0.0)


def rgb_image(img_b, img_v, percentile=99.5, stretch=0.5, q=8.0, filename=None):
    """Imagen RGB (uint8, N x M x 3) con el estiramiento asinh de Lupton et al. (2004)."""
    b = _sky_subtract(img_b)
    v = _sky_subtract(img_v)
    b /= np.percentile(b[b > 0], percentile)
    v /= np.percentile(v[v > 0], percentile)
    g = 0.5 * (b + v)
    return make_lupton_rgb(v, g, b, minimum=0.0, stretch=stretch, Q=q, filename=filename)


def color_map(img_b, img_v, zp_eff_b, zp_eff_v, smooth_px=5.0, min_snr=3.0):
    """Mapa de color superficial B - V (mag) suavizado.

    ``zp_eff`` = ZP - k X de cada filtro. Los píxeles con poca señal quedan
    en NaN para no mostrar ruido.
    """
    b = ndimage.gaussian_filter(_sky_subtract(img_b), smooth_px)
    v = ndimage.gaussian_filter(_sky_subtract(img_v), smooth_px)
    noise_b = np.nanstd(b[b < np.percentile(b, 50)])
    noise_v = np.nanstd(v[v < np.percentile(v, 50)])
    ok = (b > min_snr * noise_b) & (v > min_snr * noise_v)
    with np.errstate(invalid="ignore", divide="ignore"):
        bv = (zp_eff_b - 2.5 * np.log10(b)) - (zp_eff_v - 2.5 * np.log10(v))
    bv[~ok] = np.nan
    return bv
