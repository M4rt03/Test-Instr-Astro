"""Reducción de la noche del 15/10/2026 (NGC 7793) con el método de 2025 (fuera de la tarea).

Replica los pasos de `Reducción3.333.ipynb` y `reduccion y det. de objetos.ipynb` (grupo de 2025),
cambiando solo los nombres de archivo y lo que los datos de 2026 obligan a cambiar:

  1. Master bias = promedio de los bias.            -> no hay bias en 2026: se usan los darks de
                                                        120 s como "bias" (nivel cero).
  2. Corriente oscura = promedio(dark - bias) / t.   -> darks de 240 s (en 2025: 100 s).
  3. Master dark de t_ciencia = bias + corriente oscura * t.
  4. Flats: flat - bias - corriente oscura * t_flat, promedio, normalizado por el PROMEDIO.
  5. Ciencia: (cruda - master dark) / flat, y PROMEDIO de las imágenes SIN ALINEAR.
  6. * BPM.                                          -> el BPM.fit de 2025 no está: se omite.
  7. FITS con PIXSCALE = 0.36 y WCS desde OBJCTRA/OBJCTDEC (CD = ±0,0001°).
  8. Detección: detect_threshold(nsigma=1), detect_sources(npixels=10),
     deblend_sources(npixels=10, nlevels=32, contrast=1), SourceCatalog.

No se corrigen los errores del método (sin alineamiento, escala 0,36″/px, WCS del header, normalización
por el promedio, sin exclusión de imágenes malas): el objetivo es ver cómo se ven los datos del jueves
procesados así, para compararlos con el pipeline actual.

    python reduccion_2025_15oct.py
"""

import glob
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.visualization import ZScaleInterval
from photutils.segmentation import SourceCatalog, deblend_sources, detect_sources, detect_threshold

warnings.simplefilter("ignore")

AQUI = Path(__file__).parent
DATOS = AQUI / ".." / ".." / "pipeline" / "datos" / "2026"
CIENCIA = DATOS / "grupo4" / "octubre 15 2026"
PIPELINE = AQUI / ".." / ".." / "pipeline" / "resultados" / "noche_15oct"
SALIDA = AQUI / "salida"
SALIDA.mkdir(exist_ok=True)


def leer(f):
    return fits.getdata(f).astype(np.float32)


def lista(patron):
    return sorted(glob.glob(str(patron)))


# 1. "Master bias": promedio de los darks de 120 s (no hay bias)
bias_files = lista(DATOS / "darks" / "Dark_120.000secs_*.fit")
master_bias = np.mean([leer(f) for f in bias_files], axis=0)

# 2. Corriente oscura con los darks de 240 s
dark_files = lista(DATOS / "darks" / "Dark_240.000secs_*.fit")
t_dark = fits.getheader(dark_files[0])["EXPTIME"]
dark_current = np.mean([leer(f) - master_bias for f in dark_files], axis=0) / t_dark

# 4. Flats (corregidos con el tiempo de cada uno, promedio, normalizados por el promedio)
flats = {}
for filt in ("B", "V"):
    corr = [leer(f) - master_bias - dark_current * fits.getheader(f)["EXPTIME"]
            for f in lista(DATOS / "flats" / f"FLAT_{filt}_*.fit")]
    mf = np.mean(corr, axis=0)
    flats[filt] = mf / np.mean(mf)
    del corr

# 3, 5 y 7. Ciencia: todas las de 120 s, sin alinear y sin excluir; FITS como en 2025
finales, headers, archivos = {}, {}, {}
for filt in ("B", "V"):
    files = lista(CIENCIA / f"NGC_7793_{filt}_120.000secs_*.fit")
    t = fits.getheader(files[0])["EXPTIME"]
    master_dark = master_bias + dark_current * t
    red = [(leer(f) - master_dark) / flats[filt] for f in files]
    finales[filt] = np.mean(red, axis=0)
    del red
    h = fits.getheader(files[0])
    h["PIXSCALE"] = 0.36
    h["CTYPE1"], h["CTYPE2"] = "RA---TAN", "DEC--TAN"
    h["CRPIX1"], h["CRPIX2"] = int(h["NAXIS1"] / 2), int(h["NAXIS2"] / 2)
    coord = SkyCoord(f"{h['OBJCTRA']} {h['OBJCTDEC']}", unit=(u.hourangle, u.deg))
    h["CRVAL1"], h["CRVAL2"] = coord.ra.value, coord.dec.value
    h["CD1_1"], h["CD1_2"], h["CD2_1"], h["CD2_2"] = -0.0001, 0, 0, 0.0001
    fits.PrimaryHDU(header=h, data=np.float32(finales[filt])).writeto(SALIDA / f"Img_{filt}.fits", overwrite=True)
    headers[filt], archivos[filt] = h, [Path(f).name for f in files]
    print(f"{filt}: {len(files)} imágenes de {t:g} s promediadas sin alinear -> salida/Img_{filt}.fits")

# 8. Detección como en 2025
for filt in ("V", "B"):
    data = finales[filt]
    thr = detect_threshold(data, nsigma=1)
    sour = detect_sources(data, thr, npixels=10)
    deb = deblend_sources(data, sour, npixels=10, nlevels=32, contrast=1, progress_bar=False)
    tab = SourceCatalog(data, deb).to_table()["label", "xcentroid", "ycentroid", "area", "segment_flux"]
    tab.write(SALIDA / f"catalogo_2025_{filt}.ecsv", overwrite=True)
    print(f"{filt}: detect_sources {sour.nlabels} fuentes; tras deblend {deb.nlabels}; "
          f"área mediana {np.median(np.asarray(tab['area'])):.0f} px, máxima {np.max(np.asarray(tab['area'])):.0f} px")

# Figuras: como en el notebook (zscale, gris) y comparación con el pipeline actual
zs = ZScaleInterval(contrast=0.25)
fig, ax = plt.subplots(1, 2, figsize=(12, 6))
for a, filt in zip(ax, ("B", "V")):
    lo, hi = zs.get_limits(finales[filt])
    im = a.imshow(finales[filt], vmin=lo, vmax=hi, cmap="gray", origin="lower")
    a.set_title(f"Imagen reducida filtro {filt} (método 2025)")
    fig.colorbar(im, ax=a, label="Intensidad [ADU]", fraction=0.046, pad=0.04)
plt.tight_layout()
plt.savefig(AQUI / "reducidas_metodo2025_15oct.png", dpi=110)
plt.close()

# Recorte alrededor de una estrella brillante: método 2025 frente al pipeline actual
x0, y0, r = 996, 3865, 120          # estrella de V ≈ 11,4 (APASS), lejos de la galaxia
fig, ax = plt.subplots(2, 2, figsize=(9, 9))
for j, filt in enumerate(("B", "V")):
    c25 = finales[filt][y0 - r:y0 + r, x0 - r:x0 + r]
    pf = PIPELINE / f"stack_{filt}.fits"
    cpi = fits.getdata(pf)[y0 - r:y0 + r, x0 - r:x0 + r] if pf.exists() else np.zeros_like(c25)
    for i, (c, nombre) in enumerate(((c25, "método 2025 (sin alinear)"), (cpi, "pipeline actual (alineado)"))):
        lo, hi = np.nanpercentile(c, [5, 99.7])
        ax[j, i].imshow(c, vmin=lo, vmax=hi, cmap="gray", origin="lower")
        ax[j, i].set_title(f"{filt}: {nombre}", fontsize=10)
        ax[j, i].set_xticks([]), ax[j, i].set_yticks([])
fig.suptitle(f"Estrella en ({x0}, {y0}): recorte de {2 * r} × {2 * r} px ({2 * r * 0.316 / 60:.1f}′)", fontsize=11)
plt.tight_layout()
plt.savefig(AQUI / "comparacion_estrella_15oct.png", dpi=110)
plt.close()

# Dónde dice el header que está la imagen, frente a dónde está realmente NGC 7793
real = SkyCoord("23h57m49.75s -32d35m27.7s")
hdr = SkyCoord(headers["V"]["CRVAL1"], headers["V"]["CRVAL2"], unit="deg")
print(f"WCS de 2025 (desde OBJCTRA): centro a {hdr.separation(real).arcmin:.0f}′ de NGC 7793; "
      f"escala 0,36″/px (la real es 0,316″/px)")
print("Figuras: reducidas_metodo2025_15oct.png y comparacion_estrella_15oct.png")
