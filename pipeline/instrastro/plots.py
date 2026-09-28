"""Figuras para el informe."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.visualization import ZScaleInterval  # noqa: E402


def show(ax, data, title="", cmap="gray", label="ADU/s", zscale=True, fig=None):
    """Muestra una imagen con escala ZScale y barra de color."""
    finite = data[np.isfinite(data)]
    vmin, vmax = ZScaleInterval(contrast=0.25).get_limits(finite) if zscale and finite.size else (None, None)
    im = ax.imshow(data, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_title(title)
    if fig is not None:
        fig.colorbar(im, ax=ax, label=label, fraction=0.046, pad=0.04)
    return im


def grid(images, titles, path, label="ADU/s", ncols=None, size=5):
    """Guarda una fila (o grilla) de imágenes con sus títulos."""
    n = len(images)
    ncols = ncols or n
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(size * ncols, size * nrows), squeeze=False)
    for ax, img, title in zip(axes.flat, images, titles):
        show(ax, img, title, label=label, fig=fig)
    for ax in axes.flat[n:]:
        ax.axis("off")
    save(fig, path)


def segmentation(ax, segm, title):
    ax.imshow(segm.data, origin="lower", cmap=segm.make_cmap(seed=123), interpolation="nearest")
    ax.set_title(f"{title} ({len(segm.labels)} fuentes)")


def save(fig, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
