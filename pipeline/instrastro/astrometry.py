"""Calibración astrométrica.

Tres modos (clave ``astrometria.modo`` del YAML):

- ``header``: WCS aproximado a partir de las coordenadas del header
  (OBJCTRA/OBJCTDEC) y la escala de placa. Sirve para una primera mirada, pero
  no corrige la rotación de la cámara ni el error de apuntado.
- ``archivo``: se sube el apilado a https://nova.astrometry.net, se descarga
  el ``wcs.fits`` resultante y se indica su ruta en ``astrometria.wcs.<filtro>``.
- ``astrometry_net``: resuelve automáticamente con la API de astrometry.net
  (requiere ``pip install astroquery`` y una API key gratuita de nova.astrometry.net).
"""

import warnings

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import WCS

WCS_KEYS = ("WCSAXES", "CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2",
            "CD1_1", "CD1_2", "CD2_1", "CD2_2", "CDELT1", "CDELT2", "PC1_1", "PC1_2", "PC2_1", "PC2_2",
            "EQUINOX", "RADESYS", "LONPOLE", "LATPOLE", "A_ORDER", "B_ORDER", "AP_ORDER", "BP_ORDER")


def target_coord(header=None, ra=None, dec=None):
    """Coordenadas del objeto: de la configuración si se dan, si no del header."""
    if ra is None or dec is None:
        ra, dec = header["OBJCTRA"], header["OBJCTDEC"]
    if isinstance(ra, str) and (":" in ra or " " in ra.strip()):
        return SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    return SkyCoord(float(ra), float(dec), unit=u.deg)


def initial_wcs(shape, coord, pixscale_arcsec):
    """WCS TAN aproximado: norte arriba, este a la izquierda, centro en ``coord``."""
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [coord.ra.deg, coord.dec.deg]
    # CRPIX usa la convención FITS (el primer píxel es 1)
    w.wcs.crpix = [shape[1] / 2 + 0.5, shape[0] / 2 + 0.5]
    scale = pixscale_arcsec / 3600.0
    # RA crece hacia el este (izquierda): CD1_1 negativo; Dec crece hacia arriba: CD2_2 positivo
    w.wcs.cd = [[-scale, 0.0], [0.0, scale]]
    return w


def strip_wcs(header):
    """Borra las claves WCS de un header (antes de escribir una solución nueva)."""
    for key in list(header.keys()):
        if key in WCS_KEYS or key.startswith(("A_", "B_", "AP_", "BP_", "PV1_", "PV2_")):
            del header[key]
    return header


def apply_wcs(header, wcs):
    """Copia una solución WCS (objeto WCS o header) al header de la imagen."""
    if not isinstance(wcs, WCS):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", AstropyWarning)
            wcs = WCS(wcs)
    strip_wcs(header)
    header.update(wcs.to_header(relax=True))
    return header


def load_wcs_file(path):
    """Lee el ``wcs.fits`` (o ``new-image.fits``) descargado de astrometry.net."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        return WCS(fits.getheader(path))


def solve_astrometry_net(x, y, flux, shape, api_key, coord=None, pixscale_arcsec=None,
                         radius_deg=1.0, timeout=600):
    """Resuelve con la API de astrometry.net subiendo solo la lista de fuentes.

    Subir la lista (x, y) de las ~200 fuentes más brillantes es mucho más
    rápido que subir una imagen de 4096x4096.
    """
    try:
        from astroquery.astrometry_net import AstrometryNet
    except ImportError as exc:  # pragma: no cover - dependencia opcional
        raise ImportError("Instala astroquery (pip install astroquery) para usar el modo astrometry_net") from exc

    order = np.argsort(np.asarray(flux))[::-1][:200]
    settings = {}
    if coord is not None:
        settings.update(center_ra=coord.ra.deg, center_dec=coord.dec.deg, radius=radius_deg)
    if pixscale_arcsec is not None:
        settings.update(scale_units="arcsecperpix", scale_type="ev", scale_est=pixscale_arcsec, scale_err=10)
    ast = AstrometryNet()
    ast.api_key = api_key
    # astrometry.net usa coordenadas de píxel FITS (empiezan en 1)
    header = ast.solve_from_source_list(np.asarray(x)[order] + 1, np.asarray(y)[order] + 1,
                                        shape[1], shape[0], solve_timeout=timeout, **settings)
    if not header:
        raise RuntimeError("astrometry.net no encontró solución")
    return WCS(header)


def wcs_summary(wcs, shape):
    """Escala de placa (arcsec/px), rotación (grados) y si la imagen está espejada."""
    cd = wcs.pixel_scale_matrix
    det = np.linalg.det(cd)
    scale = np.sqrt(abs(det)) * 3600
    rotation = np.degrees(np.arctan2(cd[1, 0], cd[1, 1]))
    center = wcs.pixel_to_world_values((shape[1] - 1) / 2.0, (shape[0] - 1) / 2.0)
    return {
        "escala_arcsec_px": float(scale),
        "rotacion_deg": float(rotation),
        "espejada": bool(det > 0),
        "ra_centro_deg": float(center[0]),
        "dec_centro_deg": float(center[1]),
        "CD": cd.tolist(),
    }
