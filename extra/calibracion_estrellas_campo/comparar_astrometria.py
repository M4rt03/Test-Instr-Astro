"""Compara la astrometría de astrometry.net (header del apilado) con la ajustada contra Gaia.

    python comparar_astrometria.py ../../pipeline/resultados/noche_15oct
"""
import sys
import warnings

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

import calibrar_campo as cc

warnings.simplefilter("ignore")
carpeta = sys.argv[1]
data, _, hdr = cc.leer_stack(carpeta, "V")
w_anet = WCS(hdr)
ny, nx = data.shape
obj = SkyCoord.from_name(str(hdr.get("OBJECT", "")).strip())
det = cc.detectar(data, fwhm_px=12.0)[:400]
x, y = np.asarray(det["xcentroid"]), np.asarray(det["ycentroid"])
cat = cc.catalogo_gaia(obj, cc.ESCALA * max(nx, ny) / 60 / np.sqrt(2) + 2)
w_gaia, idx, ok = cc.resolver_astrometria(w_anet, x, y, cat, obj)
cg = SkyCoord(cat["ra"], cat["dec"], unit="deg")


def resumen(w):
    cd = w.pixel_scale_matrix
    escala = np.sqrt(abs(np.linalg.det(cd))) * 3600
    rot = np.degrees(np.arctan2(cd[1, 0], cd[1, 1]))
    return escala, rot, w.pixel_to_world(nx / 2, ny / 2)


print("\n                      escala (\"/px)   rotación (°)   centro de la imagen")
for nombre, w in (("astrometry.net", w_anet), ("ajuste con Gaia", w_gaia)):
    e, r, c = resumen(w)
    print(f"{nombre:18s}   {e:.5f}         {r:+.3f}        {c.to_string('hmsdms', precision=2)}")

# Diferencia entre ambas soluciones en una grilla de la imagen
gx, gy = np.meshgrid(np.linspace(0, nx - 1, 9), np.linspace(0, ny - 1, 9))
d = w_anet.pixel_to_world(gx.ravel(), gy.ravel()).separation(w_gaia.pixel_to_world(gx.ravel(), gy.ravel())).arcsec
print(f"\nDiferencia entre las dos soluciones en la imagen: mediana {np.median(d):.2f}\", máxima {d.max():.2f}\"")

# Residuos de cada solución contra Gaia, con las mismas estrellas
for nombre, w in (("astrometry.net", w_anet), ("ajuste con Gaia", w_gaia)):
    s = w.pixel_to_world(x[ok], y[ok])
    r = s.separation(cg[idx[ok]]).arcsec
    print(f"Residuo contra Gaia ({nombre}): rms {np.sqrt(np.mean(r ** 2)):.2f}\" con {ok.sum()} estrellas")
print(f"\nCentro de NGC 7793 (SIMBAD) en píxeles según astrometry.net: "
      f"{np.round(w_anet.world_to_pixel(obj), 0)}")
