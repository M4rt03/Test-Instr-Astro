"""Lectura de la configuración y de archivos FITS."""

import fnmatch
import glob
import warnings
from pathlib import Path

import numpy as np
import yaml
from astropy import units as u
from astropy.io import fits
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning


def load_config(path):
    """Lee el archivo YAML de configuración.

    Las rutas relativas dentro del YAML se interpretan respecto de la carpeta
    donde está el propio YAML (clave interna ``_base``).
    """
    path = Path(path)
    with open(path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    cfg["_base"] = path.parent.resolve()
    return cfg


def resolve(cfg, rel):
    """Convierte una ruta del YAML en ruta absoluta."""
    p = Path(rel)
    return p if p.is_absolute() else Path(cfg["_base"]) / p


def find_files(cfg, patterns, exclude=()):
    """Busca archivos con uno o varios patrones glob, descartando los de ``exclude``.

    ``exclude`` son patrones sobre el nombre del archivo, por ejemplo
    ``"*00000529*"`` para sacar una imagen con un satélite.
    """
    if not patterns:
        return []
    if isinstance(patterns, str):
        patterns = [patterns]
    found = set()
    for pat in patterns:
        found.update(glob.glob(str(resolve(cfg, pat))))
    keep = [f for f in found if not any(fnmatch.fnmatch(Path(f).name, ex) for ex in exclude or ())]
    return [Path(f) for f in sorted(keep)]


def read_fits(path):
    """Lee la primera HDU con datos. Devuelve ``(data float32, header)``."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        with fits.open(path) as hdul:
            hdu = next(h for h in hdul if h.data is not None)
            return np.asarray(hdu.data, dtype=np.float32), hdu.header.copy()


def write_fits(path, data, header=None, masks=None):
    """Escribe una imagen float32 y, opcionalmente, máscaras como extensiones uint8."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    hdus = [fits.PrimaryHDU(np.asarray(data, dtype=np.float32), header=header)]
    for name, mask in (masks or {}).items():
        hdus.append(fits.ImageHDU(np.asarray(mask, dtype=np.uint8), name=name))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        fits.HDUList(hdus).writeto(path, overwrite=True)


def exptime(header):
    """Tiempo de exposición en segundos."""
    for key in ("EXPTIME", "EXPOSURE"):
        if key in header:
            return float(header[key])
    raise KeyError("El header no tiene EXPTIME ni EXPOSURE")


def filter_name(header):
    """Nombre del filtro según el header (o None)."""
    value = header.get("FILTER")
    return str(value).strip() if value is not None else None


def obs_time(header):
    """Instante medio de la exposición (``DATE-OBS`` + EXPTIME/2), en UTC."""
    date = str(header["DATE-OBS"])
    if "T" not in date and "TIME-OBS" in header:
        date = f"{date}T{header['TIME-OBS']}"
    t = Time(date, format="isot", scale="utc")
    try:
        return t + exptime(header) / 2 * u.s
    except KeyError:
        return t
