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


def _stars(img, n=150, sigma=2.5, box=7):
    """Posiciones ``(x, y)`` de las ``n`` estrellas más brillantes (máximos locales + centroide)."""
    a = np.array(img, dtype=np.float64)
    good = np.isfinite(a)
    a[~good] = np.median(a[good]) if good.any() else 0.0
    a = a - ndimage.gaussian_filter(a, 25)                # quita la galaxia o nebulosa
    s = ndimage.gaussian_filter(a, sigma)
    med = np.median(s[::4, ::4])
    mad = 1.4826 * np.median(np.abs(s[::4, ::4] - med))
    peaks = (s == ndimage.maximum_filter(s, size=15)) & (s > med + 10 * mad)
    ys, xs = np.nonzero(peaks)
    edge = 3 * box
    keep = (xs > edge) & (xs < a.shape[1] - edge) & (ys > edge) & (ys < a.shape[0] - edge)
    ys, xs = ys[keep], xs[keep]
    order = np.argsort(s[ys, xs])[::-1][:n]
    pos = []
    for y, x in zip(ys[order], xs[order]):
        c = np.clip(a[y - box:y + box + 1, x - box:x + box + 1], 0, None)
        if c.sum() > 0:
            cy, cx = ndimage.center_of_mass(c)
            pos.append((x - box + cx, y - box + cy))
    return np.array(pos)


def measure_shift_stars(ref, img, max_shift=300, min_matches=5):
    """Desplazamiento ``(dy, dx)`` con las posiciones de las estrellas, o ``None`` si no alcanza.

    Vota la diferencia de posición entre todas las parejas de estrellas (la galaxia no
    interviene) y refina con la mediana de las parejas que coinciden.
    """
    p, q = _stars(ref), _stars(img)
    if len(p) < min_matches or len(q) < min_matches:
        return None
    d = (p[:, None, :] - q[None, :, :]).reshape(-1, 2)
    d = d[np.all(np.abs(d) < max_shift, axis=1)]
    if len(d) == 0:
        return None
    h, xe, ye = np.histogram2d(d[:, 0], d[:, 1], bins=2 * max_shift // 3,
                               range=[[-max_shift, max_shift], [-max_shift, max_shift]])
    i, j = np.unravel_index(np.argmax(h), h.shape)
    if h[i, j] < min_matches:
        return None
    dx0, dy0 = (xe[i] + xe[i + 1]) / 2, (ye[j] + ye[j + 1]) / 2
    sel = (np.abs(d[:, 0] - dx0) < 4) & (np.abs(d[:, 1] - dy0) < 4)
    dx, dy = np.median(d[sel], axis=0)
    return float(dy), float(dx)


def measure_shift(ref, img, size=2048):
    """Desplazamiento ``(dy, dx)`` que hay que aplicar a ``img`` para alinearla con ``ref``.

    Primero usa las posiciones de las estrellas (robusto frente a una galaxia o nebulosa que
    llena el centro). Si hay muy pocas estrellas, usa la correlación cruzada (vía FFT) de la
    región central de ``size`` x ``size`` píxeles, refinada con un ajuste parabólico.
    """
    a = _prepare(ref, size)
    b = _prepare(img, size)
    corr = np.fft.ifft2(np.fft.fft2(a) * np.conj(np.fft.fft2(b))).real
    n = corr.shape[0]
    coarse = measure_shift_stars(ref, img)
    if coarse is not None:
        # Las estrellas dan el desplazamiento aproximado; la correlación lo refina a sub-píxel,
        # buscando el máximo solo a ±3 px de ese valor (así la galaxia no desvía el máximo).
        cy, cx = int(round(coarse[0])) % n, int(round(coarse[1])) % n
        win = np.full_like(corr, -np.inf)
        for oy in range(-3, 4):
            for ox in range(-3, 4):
                win[(cy + oy) % n, (cx + ox) % n] = corr[(cy + oy) % n, (cx + ox) % n]
        iy, ix = np.unravel_index(np.argmax(win), win.shape)
    else:
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
