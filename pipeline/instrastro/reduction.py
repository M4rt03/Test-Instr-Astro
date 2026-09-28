"""Reducción de imágenes: master bias, dark current, master flat y apilado.

Convenciones:

- Todas las imágenes reducidas quedan en **ADU/s** (se dividen por su propio
  tiempo de exposición). Así se pueden combinar exposiciones de distinta
  duración y la fotometría no necesita volver a dividir por EXPTIME.
- Los píxeles malos (BPM, esquinas viñeteadas, bordes tras alinear) se marcan
  con ``NaN``, no con 0: un 0 contamina las medianas, el fondo y las
  magnitudes; un NaN simplemente se ignora.
- Los píxeles saturados en la imagen cruda se guardan en una máscara aparte,
  para excluir esas estrellas de la fotometría y del seeing.
"""

import warnings
from dataclasses import dataclass, field

import numpy as np
from astropy.stats import sigma_clip
from astropy.table import Table

from . import alignment
from .io import exptime, filter_name, read_fits


def _central(img, frac=0.5):
    ny, nx = img.shape
    dy, dx = int(ny * frac / 2), int(nx * frac / 2)
    return img[ny // 2 - dy:ny // 2 + dy, nx // 2 - dx:nx // 2 + dx]


def combine(frames, method="mean", sigma=3.0, chunk_rows=256):
    """Combina una lista de imágenes píxel a píxel ignorando NaN.

    ``method="mean"`` hace un promedio con rechazo sigma-clipping (elimina
    rayos cósmicos y satélites y mantiene el S/N del promedio). Con menos de 4
    imágenes el clipping no es confiable y se usa la mediana. Se procesa por
    bloques de filas para no ocupar toda la RAM con imágenes de 4096x4096.
    """
    frames = list(frames)
    if not frames:
        raise ValueError("No hay imágenes para combinar")
    if len(frames) == 1:
        return np.array(frames[0], dtype=np.float32)
    ny, nx = frames[0].shape
    out = np.empty((ny, nx), dtype=np.float32)
    use_median = method == "median" or len(frames) < 4
    for r0 in range(0, ny, chunk_rows):
        cube = np.stack([f[r0:r0 + chunk_rows] for f in frames]).astype(np.float32)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # NaN en bordes y píxeles malos: se ignoran a propósito
            if use_median:
                out[r0:r0 + chunk_rows] = np.nanmedian(cube, axis=0)
            else:
                clipped = sigma_clip(cube, sigma=sigma, maxiters=3, axis=0, cenfunc="median",
                                     stdfunc="mad_std", masked=False, copy=False)
                out[r0:r0 + chunk_rows] = np.nanmean(clipped, axis=0)
    return out


def make_master_bias(files):
    """Master bias = mediana de los bias (robusta a rayos cósmicos)."""
    return combine([read_fits(f)[0] for f in files], method="median")


def make_dark_current(files, master_bias):
    """Corriente oscura en ADU/s/píxel: promedio de (dark - bias) / t_dark.

    Al separar la parte que depende del tiempo (dark current) de la que no
    (bias), el dark se puede escalar a cualquier tiempo de exposición.
    """
    frames = []
    for f in files:
        data, hdr = read_fits(f)
        frames.append((data - master_bias) / exptime(hdr))
    return combine(frames, method="mean")


def make_master_flat(files, master_bias, dark_current, saturation=None, min_level=0.2):
    """Master flat normalizado (mediana = 1 en la zona central).

    Cada flat se corrige por bias + dark escalado a su propio tiempo y se
    normaliza por su propia mediana **antes** de combinar: los flats de
    crepúsculo cambian de brillo de una toma a otra, y promediar flats de
    distinto nivel sin normalizar pesa más a los más brillantes.
    Los píxeles con respuesta < ``min_level`` (esquinas viñeteadas) quedan en NaN.

    Devuelve ``(flat, niveles)``, donde ``niveles`` es la mediana de cada flat
    en ADU (útil para comprobar que no estén saturados ni subexpuestos).
    """
    frames, levels = [], []
    for f in files:
        data, hdr = read_fits(f)
        corr = data - master_bias - dark_current * exptime(hdr)
        level = float(np.nanmedian(_central(corr)))
        levels.append(level)
        if saturation and level > 0.8 * saturation:
            warnings.warn(f"{f.name}: nivel {level:.0f} ADU cerca de saturación; no es lineal")
        frames.append(corr / level)
    flat = combine(frames, method="median")
    flat /= np.nanmedian(_central(flat))
    flat[~(flat > min_level)] = np.nan
    return flat, levels


@dataclass
class Calibration:
    """Imágenes maestras de una noche."""

    bias: np.ndarray
    dark_current: np.ndarray
    flats: dict
    bad_pixels: np.ndarray = None
    saturation: float = 65000.0
    flat_levels: dict = field(default_factory=dict)


def reduce_frame(raw, header, calib, filt=None):
    """Reduce una imagen cruda. Devuelve ``(reducida en ADU/s, máscara de saturación)``.

    reducida = (cruda - bias - dark_current * t) / flat / t
    """
    filt = filt or filter_name(header)
    if filt not in calib.flats:
        raise KeyError(f"No hay master flat para el filtro {filt!r}")
    t = exptime(header)
    saturated = raw >= calib.saturation
    with np.errstate(invalid="ignore", divide="ignore"):
        red = (raw - calib.bias - calib.dark_current * t) / calib.flats[filt] / t
    if calib.bad_pixels is not None:
        red[calib.bad_pixels] = np.nan
    return red.astype(np.float32), saturated


def reduce_and_stack(files, calib, filt, reference=None, sigma=3.0, measure=None, airmass_func=None):
    """Reduce, alinea y combina una lista de imágenes de ciencia de un filtro.

    Parámetros
    ----------
    reference : imagen a la que se alinean todas (por defecto, la primera reducida).
    measure : función opcional ``f(reducida, saturadas) -> dict`` que se aplica a
        cada exposición individual antes de alinear (se usa para medir el seeing
        de cada exposición).
    airmass_func : función opcional ``f(header) -> masa de aire``.

    Devuelve ``(apilado, máscara_saturación, header, tabla_por_exposición)``.
    """
    reduced, saturated, headers, rows = [], [], [], []
    for f in files:
        raw, hdr = read_fits(f)
        if filter_name(hdr) and filter_name(hdr) != filt:
            warnings.warn(f"{f.name}: el header dice filtro {filter_name(hdr)!r}, se esperaba {filt!r}")
        red, sat = reduce_frame(raw, hdr, calib, filt)
        row = {"archivo": f.name, "exptime": exptime(hdr), "fondo_adu_s": float(np.nanmedian(red))}
        if airmass_func is not None:
            row["airmass"] = airmass_func(hdr)
        if measure is not None:
            row.update(measure(red, sat))
        reduced.append(red)
        saturated.append(sat)
        headers.append(hdr)
        rows.append(row)

    ref = reduced[0] if reference is None else reference
    for i, red in enumerate(reduced):
        dy, dx = alignment.measure_shift(ref, red)
        rows[i]["shift_y"], rows[i]["shift_x"] = dy, dx
        reduced[i] = alignment.shift_image(red, dy, dx)
        saturated[i] = alignment.shift_mask(saturated[i], dy, dx)

    # La interpolación bilineal promedia píxeles vecinos: baja la varianza por píxel en un factor
    # ((1-fx)^2 + fx^2)((1-fy)^2 + fy^2) pero no la del flujo sumado en una apertura. NOISECOR
    # corrige la desviación estándar medida por píxel en el apilado para usarla en errores.
    reduction_factor = []
    for row in rows:
        fy, fx = abs(row["shift_y"]) % 1, abs(row["shift_x"]) % 1
        reduction_factor.append(((1 - fx) ** 2 + fx ** 2) * ((1 - fy) ** 2 + fy ** 2))
    noise_corr = float(np.sqrt(len(reduction_factor) / np.sum(reduction_factor)))

    stack = combine(reduced, method="mean", sigma=sigma)
    sat_mask = np.logical_or.reduce(saturated)
    info = Table(rows=rows)

    header = headers[0].copy()
    total = float(np.sum(info["exptime"]))
    header["BUNIT"] = ("ADU/s", "imagen reducida y normalizada por tiempo")
    header["NCOMBINE"] = (len(files), "imagenes combinadas")
    header["TOTEXP"] = (total, "[s] tiempo total de exposicion")
    header["FILTER"] = filt
    header["NOISECOR"] = (noise_corr, "factor ruido correlacionado por el alineamiento")
    if "airmass" in info.colnames:
        header["AIRMASS"] = (float(np.nanmean(info["airmass"])), "promedio de las exposiciones")
    return stack, sat_mask, header, info
