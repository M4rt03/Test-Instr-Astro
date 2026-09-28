"""Alineamiento de imágenes por correlación cruzada.

El telescopio deriva unos pocos píxeles entre exposiciones (y entre filtros),
así que antes de combinar hay que llevar todas las imágenes a la misma grilla.
Si no se alinean, las estrellas del apilado quedan "corridas" y el seeing
medido sale artificialmente peor.
"""

import numpy as np
from scipy import ndimage


def _prepare(img, size):
    """Recorta el centro, quita el fondo extendido (nebulosa) y limita los picos."""
    ny, nx = img.shape
    s = min(size, ny, nx)
    y0, x0 = (ny - s) // 2, (nx - s) // 2
    a = np.array(img[y0:y0 + s, x0:x0 + s], dtype=np.float64)
    good = np.isfinite(a)
    a[~good] = np.median(a[good]) if good.any() else 0.0
    # Filtro pasa-altos: deja solo las estrellas, que son las que fijan la posición
    a -= ndimage.gaussian_filter(a, 15)
    a = np.clip(a, 0, np.percentile(a, 99.9))
    window = np.outer(np.hanning(s), np.hanning(s))
    return a * window


def _parabola(cm, c0, cp):
    """Vértice de la parábola que pasa por tres puntos equiespaciados."""
    denom = cm - 2 * c0 + cp
    return 0.0 if denom == 0 else 0.5 * (cm - cp) / denom


def measure_shift(ref, img, size=2048):
    """Desplazamiento ``(dy, dx)`` que hay que aplicar a ``img`` para alinearla con ``ref``.

    Usa la correlación cruzada (vía FFT) de la región central de ``size`` x ``size``
    píxeles y refina el máximo a precisión sub-píxel con un ajuste parabólico.
    """
    a = _prepare(ref, size)
    b = _prepare(img, size)
    corr = np.fft.ifft2(np.fft.fft2(a) * np.conj(np.fft.fft2(b))).real
    n = corr.shape[0]
    iy, ix = np.unravel_index(np.argmax(corr), corr.shape)
    fy = _parabola(corr[(iy - 1) % n, ix], corr[iy, ix], corr[(iy + 1) % n, ix])
    fx = _parabola(corr[iy, (ix - 1) % n], corr[iy, ix], corr[iy, (ix + 1) % n])
    dy = iy + fy
    dx = ix + fx
    # La correlación es periódica: desplazamientos mayores a n/2 son negativos
    if dy > n / 2:
        dy -= n
    if dx > n / 2:
        dx -= n
    return float(dy), float(dx)


def shift_image(img, dy, dx):
    """Desplaza la imagen con interpolación bilineal; lo que queda fuera se marca como NaN."""
    bad = ~np.isfinite(img)
    filled = np.where(bad, 0.0, img)
    out = ndimage.shift(filled, (dy, dx), order=1, mode="constant", cval=np.nan)
    bad_shifted = ndimage.shift(bad.astype(np.float32), (dy, dx), order=1, mode="constant", cval=1.0)
    out[bad_shifted > 0.01] = np.nan
    return out.astype(np.float32)


def shift_mask(mask, dy, dx):
    """Desplaza una máscara booleana (un píxel queda marcado si toca uno marcado)."""
    shifted = ndimage.shift(mask.astype(np.float32), (dy, dx), order=1, mode="constant", cval=0.0)
    return shifted > 0.01
