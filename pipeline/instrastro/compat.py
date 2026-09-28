"""Compatibilidad entre versiones de photutils.

photutils 3.0 renombró varios argumentos (``npixels`` -> ``n_pixels``,
``nlevels`` -> ``n_levels``, ...) y columnas del catálogo (``xcentroid`` ->
``x_centroid``). Estas funciones detectan la versión instalada para que el
pipeline funcione tanto con las versiones antiguas (1.x, 2.x) como con la 3.x.
"""

import inspect

from photutils import segmentation as _seg


def _params(func):
    return inspect.signature(func).parameters


def _pick(params, *names):
    for name in names:
        if name in params:
            return name
    return None


def detect_sources(data, threshold, n_pixels, mask=None):
    """Envoltorio de ``photutils.segmentation.detect_sources``."""
    params = _params(_seg.detect_sources)
    kwargs = {_pick(params, "n_pixels", "npixels"): n_pixels}
    if mask is not None:
        kwargs["mask"] = mask
    return _seg.detect_sources(data, threshold, **kwargs)


def deblend_sources(data, segm, n_pixels, n_levels=32, contrast=0.001):
    """Envoltorio de ``photutils.segmentation.deblend_sources`` (sin barra de progreso)."""
    params = _params(_seg.deblend_sources)
    kwargs = {
        _pick(params, "n_pixels", "npixels"): n_pixels,
        _pick(params, "n_levels", "nlevels"): n_levels,
        "contrast": contrast,
    }
    if "progress_bar" in params:
        kwargs["progress_bar"] = False
    return _seg.deblend_sources(data, segm, **kwargs)


def cat_value(cat, *names):
    """Devuelve la primera propiedad existente del ``SourceCatalog`` entre ``names``."""
    for name in names:
        if hasattr(cat, name):
            value = getattr(cat, name)
            return getattr(value, "value", value)
    raise AttributeError(f"El catálogo no tiene ninguna de las propiedades {names}")
