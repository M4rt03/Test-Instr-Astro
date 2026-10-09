"""Revisión rápida de una exposición de prueba en el telescopio.

Para cada imagen muestra:
- header: objeto, filtro, t_exp, hora (UTC y local), temperatura del sensor, ganancia y sitio;
- masa de aire calculada para El Sauce (no la del header, que en 2023 venía con otro sitio);
- fondo del cielo en ADU, ADU/s y mag/arcsec², usando el ZP de 2025;
- máximo del centro (núcleo del objeto) y píxeles saturados, con el t_exp máximo sin saturar;
- seeing (FWHM) y elipticidad de las estrellas, para revisar foco y guiado;
- si el objeto es una estrella estándar conocida: su flujo, su pico y el ZP de la prueba.

No modifica nada. La imagen no se corrige por flat, así que el ZP y el fondo son aproximados
(error de unos pocos %).

Uso (desde la carpeta pipeline):

    python revisar_prueba.py datos/noche1/prueba_V.fit
    python revisar_prueba.py datos/noche1/*.fit
    python revisar_prueba.py prueba.fit --estandar "HD 8130"     # si OBJECT no dice el nombre
    python revisar_prueba.py prueba.fit --bias datos/noche1/Bias_1x1_00000001.fit

Para probarlo con los datos de 2025 (deberían salir ZP de ~21,97–21,99 en B para HIP 116375,
~21,70 para HIP 117678 y ~21,94 en V para las dos):

    python revisar_prueba.py datos/2025-10-02/**/V462_B_5.000secs_00000557.fit
"""

import argparse
import glob
import math
import sys
import warnings
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation
from astropy.stats import sigma_clipped_stats
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))

from instrastro import astrometry, photometry, seeing  # noqa: E402
from instrastro.io import exptime, filter_name, obs_time, read_fits  # noqa: E402

SITIO = {"lat": -30.4597, "lon": -70.7503, "alt": 1600}    # El Sauce (header de 2025)
ESCALA = 0.317           # arcsec/px (medida con Gaia)
GANANCIA = 0.85          # e-/ADU
BIAS_2025 = 92.5         # ADU, mediana del master bias de 2025
LIMITE = 40000           # ADU sobre el bias: límite de trabajo (lineal y con margen)
SATURA = 65000           # ADU crudos: saturación del convertidor
ZP_2025 = {"B": 21.845, "V": 21.941}     # ADU/s sobre la atmósfera (zeropoints.ecsv de 2025)
K = {"B": 0.25, "V": 0.15}               # extinción, mag por masa de aire
CIELO_OSCURO = {"B": 22.7, "V": 21.8}    # mag/arcsec2 en el cenit sin Luna

# Magnitudes Johnson (B, V): estrellas_estandar.py (2026) y las estándares de 2025 (Tycho
# convertido, como en el config de 2025, para poder probar el script con esos datos)
ESTANDARES = {
    "HD 202941": (7.072, 7.070), "HD 207480": (7.190, 7.140), "HD 210300": (6.590, 6.440),
    "HD 212643": (6.260, 6.290), "HD 220881": (7.720, 7.440), "HD 562": (7.787, 7.650),
    "HD 8130": (7.493, 7.448), "HD 12206": (6.810, 6.790), "HD 223884": (6.430, 6.240),
    "HD 225200": (6.386, 6.380), "HD 7323": (7.928, 7.830), "HD 195500": (7.380, 7.320),
    "HD 215863": (7.836, 7.690),
    "HIP 116375": (8.150, 7.077), "HIP 117678": (8.019, 7.058),
}


def _clave(nombre):
    return "".join(str(nombre).upper().split())


def buscar_estandar(nombre):
    tabla = {_clave(k): (k, v) for k, v in ESTANDARES.items()}
    return tabla.get(_clave(nombre))


def masa_de_aire(header):
    try:
        coord = astrometry.target_coord(header)
    except (KeyError, ValueError):
        return None, None
    lugar = EarthLocation(lat=SITIO["lat"] * u.deg, lon=SITIO["lon"] * u.deg, height=SITIO["alt"] * u.m)
    altaz = coord.transform_to(AltAz(obstime=obs_time(header), location=lugar))
    alt = float(altaz.alt.deg)
    return (float(altaz.secz) if alt > 5 else None), alt


def aviso(ok):
    return "OK " if ok else "OJO"


def revisar(path, args):
    raw, h = read_fits(path)
    t = exptime(h)
    filt = (filter_name(h) or "?").upper()
    print(f"\n=== {Path(path).name} ===")

    # --- Header
    fecha = str(h.get("DATE-OBS", "?"))
    local = ""
    try:
        hh, mm = int(fecha[11:13]), int(fecha[14:16])
        local = f" (local {(hh - 3) % 24:02d}:{mm:02d})"
    except ValueError:
        pass
    print(f"Objeto {h.get('OBJECT', '?')!s:<14} tipo {h.get('IMAGETYP', '?')!s:<12} "
          f"filtro {filt:<3} t_exp {t:g} s   {fecha[:19]} UTC{local}")

    temp, consigna = h.get("CCD-TEMP"), h.get("SET-TEMP")
    if temp is not None and consigna is not None:
        print(f"  [{aviso(abs(temp - consigna) <= 0.5)}] Sensor {temp:.1f} °C, consigna {consigna:.1f} °C")
    gain = h.get("GAIN") or h.get("EGAIN")
    if gain is not None:
        print(f"  [{aviso(abs(float(gain) - GANANCIA) < 0.02)}] Ganancia {float(gain):.3f} e-/ADU "
              f"(0,85 = modo 16-bit HDR)")
    lat, lon = h.get("OBSGEO-B"), h.get("OBSGEO-L")
    if lat is not None and lon is not None:
        ok = abs(float(lat) - SITIO["lat"]) < 0.1 and abs(float(lon) - SITIO["lon"]) < 0.1
        print(f"  [{aviso(ok)}] Sitio del header {float(lat):.3f}°, {float(lon):.3f}° "
              f"(El Sauce: {SITIO['lat']}°, {SITIO['lon']}°)")
    x_air, alt = masa_de_aire(h)
    if x_air is not None:
        print(f"        Altura {alt:.0f}°, masa de aire {x_air:.3f} (calculada para El Sauce; "
              f"header: {h.get('AIRMASS', '—')})")

    # --- Bias
    if args.bias:
        bias_img, _ = read_fits(args.bias)
        bias_txt = f"bias de {Path(args.bias).name}"
    else:
        bias_img = BIAS_2025
        bias_txt = f"bias {BIAS_2025} ADU (2025)"
    sat = raw >= SATURA
    data = raw - bias_img

    # --- Fondo
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, fondo, ruido = sigma_clipped_stats(data[::4, ::4], sigma=3.0, maxiters=5)
    print(f"  Fondo {fondo:.1f} ADU sobre el {bias_txt} ({fondo / t:.3f} ADU/s), "
          f"ruido {ruido:.1f} ADU por píxel")
    if filt in ZP_2025 and x_air is not None and fondo > 0:
        mu = ZP_2025[filt] - K[filt] * x_air - 2.5 * math.log10(fondo / t / ESCALA ** 2)
        print(f"        Cielo ≈ {mu:.1f} mag/arcsec² en {filt} (sin Luna en el cenit: "
              f"{CIELO_OSCURO[filt]}; con la Luna al 88 % se espera ~19)")
    elif fondo <= 0:
        print("  [OJO] Fondo ≤ 0: revisar el nivel de bias (¿otro modo de lectura?)")

    # --- Máximo del centro (núcleo) y saturación
    ny, nx = data.shape
    c = args.caja // 2
    centro = ndimage.median_filter(data[ny // 2 - c:ny // 2 + c, nx // 2 - c:nx // 2 + c], size=3)
    pico = float(np.nanmax(centro) - fondo)
    print(f"  [{aviso(pico < LIMITE)}] Máximo en el centro ({args.caja} px) {pico:.0f} ADU sobre el "
          f"fondo (límite {LIMITE})" + (f": t_exp máximo ≈ {t * LIMITE / pico:.0f} s" if pico > 0 else ""))
    n_sat = int(sat.sum())
    print(f"  [{aviso(n_sat == 0)}] Píxeles saturados (≥ {SATURA} ADU): {n_sat}")

    # --- Seeing y forma de las estrellas
    fwhm_guess = args.seeing / ESCALA
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tbl = seeing.measure_psf(data.astype(np.float32), sat, fwhm_guess, max_stars=60)
    res = seeing.summarize(tbl, ESCALA)
    if res["n_estrellas"] >= 3:
        ok_e = res["elipticidad"] < 0.15
        print(f"  [{aviso(res['seeing_arcsec'] < 4)}] Seeing {res['seeing_arcsec']:.2f}″ "
              f"({res['fwhm_px']:.1f} px, {res['n_estrellas']} estrellas)")
        print(f"  [{aviso(ok_e)}] Elipticidad {res['elipticidad']:.2f}"
              + ("" if ok_e else " — estrellas alargadas: revisar guiado o seguimiento"))
    else:
        print("        Pocas estrellas aisladas para medir el seeing (normal en estándares de 1–3 s)")

    # --- Estrella estándar
    est = buscar_estandar(args.estandar or h.get("OBJECT", ""))
    if est and filt in ("B", "V"):
        nombre, (mb, mv) = est
        m = mb if filt == "B" else mv
        red = data / t
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = photometry.measure_standard(red, sat, fwhm_guess, aperture_factor=3.0)
        pico_e = float(np.nanmax(data[int(r["y"]) - 3:int(r["y"]) + 4, int(r["x"]) - 3:int(r["x"]) + 4]) - fondo)
        print(f"  Estándar {nombre} ({filt} = {m:.3f}): FWHM {r['fwhm_px'] * ESCALA:.2f}″, "
              f"pico {pico_e:.0f} ADU" + (" SATURADA" if r["saturated"] else ""))
        if pico_e > 0:
            print(f"  [{aviso(pico_e < LIMITE and not r['saturated'])}] t_exp máximo sin pasar "
                  f"{LIMITE} ADU ≈ {t * LIMITE / pico_e:.1f} s (usar ~60 %: {0.6 * t * LIMITE / pico_e:.1f} s)")
        if r["flux"] > 0 and x_air is not None:
            zp = photometry.zeropoint(r["flux"], m, K[filt], x_air)
            print(f"  ZP ≈ {zp:.3f} ADU/s (2025: {ZP_2025[filt]:.3f}; diferencia {zp - ZP_2025[filt]:+.2f} mag)")
            print(f"        Para recalcular: python ../planificacion/tiempos_exposicion.py "
                  f"--zp-{filt.lower()} {zp - 0.18:.3f}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("archivos", nargs="+", help="imágenes FITS (se aceptan comodines)")
    p.add_argument("--estandar", help="nombre de la estándar si OBJECT no lo dice (ej. \"HD 8130\")")
    p.add_argument("--bias", help="un bias o master bias de la noche (si no, se usan 92,5 ADU)")
    p.add_argument("--seeing", type=float, default=2.5, help="seeing inicial en arcsec (2,5)")
    p.add_argument("--caja", type=int, default=600, help="lado de la caja central en píxeles (600 = 3,2′)")
    args = p.parse_args()

    archivos = []
    for patron in args.archivos:
        archivos += sorted(glob.glob(patron, recursive=True)) or [patron]
    for f in archivos:
        try:
            revisar(f, args)
        except Exception as e:      # una imagen mala no detiene la revisión de las demás
            print(f"\n=== {f} ===\n  ERROR: {e}")


if __name__ == "__main__":
    main()
