"""Calibración fotométrica con estrellas del campo (fuera de la tarea).

Alternativa para una noche sin estrellas estándar: se usan las estrellas del propio campo con
magnitudes B y V de catálogo para obtener el zeropoint (ZP) y el término de color.

Catálogos (con astroquery):
- APASS DR9 (VizieR II/336): B y V de Johnson; confiable entre V ≈ 10 y 16.
- Gaia DR3, fotometría sintética estandarizada (tabla gaiadr3.synthetic_photometry_gspc):
  B y V de Johnson calculados desde los espectros BP/RP. Se usa como comparación.

Entrada: los apilados del pipeline (resultados/<noche>/stack_B.fits y stack_V.fits), en ADU/s,
en la misma grilla y con una astrometría aproximada en el header.

Uso (desde esta carpeta):
    python calibrar_campo.py ../../pipeline/resultados/noche_15oct
    python calibrar_campo.py ../../pipeline/resultados/noche_15oct --catalogo gaia

Ajuste, por filtro:  m_catálogo − m_instrumental = ZP + c · (B − V)_catálogo
con m_instrumental = −2.5 log10(flujo en ADU/s). El ZP es el de la masa de aire del apilado;
también se da corregido a fuera de la atmósfera con k supuesto (no hay estándares para medirlo).
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import fit_wcs_from_points
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry
from photutils.detection import DAOStarFinder

warnings.simplefilter("ignore")

ESCALA = 0.317                     # arcsec/px
K = {"B": 0.25, "V": 0.15}         # extinción supuesta (típica de Cerro Tololo)
ZP_2025 = {"B": 21.845, "V": 21.941}   # ADU/s fuera de la atmósfera, estándares de 2025
EXCLUIR_RADIO_ARCMIN = 6.0         # alrededor del objeto (la galaxia contamina la fotometría)


def leer_stack(carpeta, filt):
    with fits.open(Path(carpeta) / f"stack_{filt}.fits") as h:
        data = h[0].data.astype(np.float32)
        hdr = h[0].header
        sat = h["SATMASK"].data.astype(bool) if "SATMASK" in h else np.zeros(data.shape, bool)
    return data, sat, hdr


def detectar(data, fwhm_px, n_sigma=8.0):
    _, med, std = sigma_clipped_stats(data[::2, ::2], sigma=3)
    # sin filtro de redondez: con viento las estrellas salen alargadas (redondez ~ -0,8)
    t = DAOStarFinder(fwhm=fwhm_px, threshold=n_sigma * std, sharplo=0.2, sharphi=0.9,
                      roundlo=-1.0, roundhi=1.0)(np.nan_to_num(data - med))
    t.sort("flux", reverse=True)
    return t


def fotometria(data, x, y, r_ap, r_in, r_out):
    pos = np.column_stack([x, y])
    ap, an = CircularAperture(pos, r=r_ap), CircularAnnulus(pos, r_in=r_in, r_out=r_out)
    ph = aperture_photometry(np.nan_to_num(data), [ap, an])
    bkg = ph["aperture_sum_1"] / an.area
    flux = np.asarray(ph["aperture_sum_0"] - bkg * ap.area)
    # error: Poisson de la fuente + ruido del fondo (en ADU/s; aproximado)
    _, _, std = sigma_clipped_stats(data[::4, ::4], sigma=3)
    err = np.sqrt(np.abs(flux) / 100.0 + ap.area * std ** 2)
    return flux, err


def catalogo_apass(centro, radio_arcmin):
    from astroquery.vizier import Vizier
    v = Vizier(columns=["RAJ2000", "DEJ2000", "Vmag", "e_Vmag", "Bmag", "e_Bmag", "nobs"], row_limit=-1)
    t = v.query_region(centro, radius=radio_arcmin * u.arcmin, catalog="II/336/apass9")[0]
    out = Table({"ra": np.asarray(t["RAJ2000"], float), "dec": np.asarray(t["DEJ2000"], float),
                 "B": np.asarray(t["Bmag"].filled(np.nan), float), "V": np.asarray(t["Vmag"].filled(np.nan), float),
                 "eB": np.asarray(t["e_Bmag"].filled(np.nan), float), "eV": np.asarray(t["e_Vmag"].filled(np.nan), float)})
    return out[np.isfinite(out["B"]) & np.isfinite(out["V"])]


def catalogo_gaia(centro, radio_arcmin):
    from astroquery.gaia import Gaia
    q = f"""
        SELECT g.ra, g.dec, s.v_jkc_mag AS V, s.b_jkc_mag AS B,
               s.v_jkc_flux_error / s.v_jkc_flux * 1.0857 AS eV,
               s.b_jkc_flux_error / s.b_jkc_flux * 1.0857 AS eB
        FROM gaiadr3.gaia_source AS g
        JOIN gaiadr3.synthetic_photometry_gspc AS s USING (source_id)
        WHERE 1 = CONTAINS(POINT('ICRS', g.ra, g.dec),
                           CIRCLE('ICRS', {centro.ra.deg}, {centro.dec.deg}, {radio_arcmin / 60}))
          AND s.v_jkc_mag < 17
    """
    t = Gaia.launch_job_async(q).get_results()
    t.rename_columns([c for c in t.colnames], [c.lower() for c in t.colnames])
    out = Table({"ra": np.asarray(t["ra"], float), "dec": np.asarray(t["dec"], float),
                 "B": np.asarray(t["b"], float), "V": np.asarray(t["v"], float),
                 "eB": np.asarray(t["eb"], float), "eV": np.asarray(t["ev"], float)})
    return out[np.isfinite(out["B"]) & np.isfinite(out["V"])]


def resolver_astrometria(wcs0, x, y, cat, c0, n_det=60, n_cat=80):
    """Resuelve la astrometría contra el catálogo: escala fija, se buscan la rotación, la
    orientación (paridad) y la traslación por votación; luego se ajusta un WCS con las coincidencias."""
    ny, nx = int(wcs0.pixel_shape[1]), int(wcs0.pixel_shape[0])
    cc = SkyCoord(cat["ra"], cat["dec"], unit="deg")
    plano = c0.skyoffset_frame()
    off = cc.transform_to(plano)
    xi, eta = off.lon.arcsec, off.lat.arcsec                       # arcsec, este y norte
    bc = np.argsort(cat["V"])[:n_cat]
    # para la búsqueda se usan estrellas lejos del centro: ahí está el objeto, y sus nudos
    # (regiones HII de una galaxia) se detectan como "estrellas" que no están en el catálogo
    lejos = np.hypot(x - nx / 2, y - ny / 2) > 700
    xs, ys = x[lejos][:n_det], y[lejos][:n_det]
    dx, dy = (xs - nx / 2) * ESCALA, (ys - ny / 2) * ESCALA
    mejor = (0, None)
    for par in (1, -1):
        for th in np.radians(np.arange(0, 360, 0.5)):
            u_ = par * dx
            px, py = np.cos(th) * u_ - np.sin(th) * dy, np.sin(th) * u_ + np.cos(th) * dy
            difx = (xi[bc][None, :] - px[:, None]).ravel()
            dify = (eta[bc][None, :] - py[:, None]).ravel()
            h, xe, ye = np.histogram2d(difx, dify, bins=200, range=[[-800, 800], [-800, 800]])
            i, j = np.unravel_index(np.argmax(h), h.shape)
            if h[i, j] > mejor[0]:
                mejor = (h[i, j], (par, th, (xe[i] + xe[i + 1]) / 2, (ye[j] + ye[j + 1]) / 2))
    votos, (par, th, tx, ty) = mejor
    print(f"Búsqueda: {int(votos)} coincidencias con paridad {par:+d}, rotación {np.degrees(th):.1f}°, "
          f"desplazamiento ({tx:+.0f}\", {ty:+.0f}\")")
    u_ = par * (x - nx / 2) * ESCALA
    px = np.cos(th) * u_ - np.sin(th) * (y - ny / 2) * ESCALA + tx
    py = np.sin(th) * u_ + np.cos(th) * (y - ny / 2) * ESCALA + ty
    sky = SkyCoord(px * u.arcsec, py * u.arcsec, frame=plano).icrs
    wcs = None
    for tol in (10.0, 5.0, 3.0):
        idx, d2d, _ = sky.match_to_catalog_sky(cc)
        ok = d2d.arcsec < tol
        if ok.sum() < 6:
            raise RuntimeError(f"Muy pocas estrellas coinciden con el catálogo ({ok.sum()})")
        wcs = fit_wcs_from_points((x[ok], y[ok]), cc[idx[ok]], projection="TAN")
        sky = wcs.pixel_to_world(x, y)
    idx, d2d, _ = sky.match_to_catalog_sky(cc)
    ok = d2d.arcsec < 3.0
    escala = np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix))) * 3600
    print(f"Astrometría: {ok.sum()} coincidencias, residuo rms {np.sqrt(np.mean(d2d.arcsec[ok] ** 2)):.2f}\", "
          f"escala {escala:.4f}\"/px")
    return wcs, idx, ok

def ajustar(dm, color, err, n_iter=5, con_color=True):
    """dm = ZP + c·color (o solo ZP si ``con_color`` es False), con recorte a 3 sigma.
    Devuelve (ZP, c, error ZP, error c, rms, máscara)."""
    ok = np.isfinite(dm) & np.isfinite(color)
    if not con_color:
        color = np.zeros_like(color)
    for _ in range(n_iter):
        A = np.column_stack([np.ones(ok.sum()), color[ok]]) if con_color else np.ones((ok.sum(), 1))
        w = 1 / np.maximum(err[ok], 0.01) ** 2
        cov = np.linalg.inv(A.T @ (A * w[:, None]))
        p = cov @ (A.T @ (w * dm[ok]))
        p = np.append(p, 0.0) if not con_color else p
        res = dm - (p[0] + p[1] * color)
        rms = 1.4826 * np.nanmedian(np.abs(res[ok] - np.nanmedian(res[ok])))
        nuevo = np.isfinite(res) & (np.abs(res) < 3 * max(rms, 0.01))
        if np.array_equal(nuevo, ok):
            break
        ok = nuevo
    escala = max(1.0, np.sqrt(np.sum((res[ok] / np.maximum(err[ok], 0.01)) ** 2) / max(ok.sum() - 2, 1)))
    ec = np.sqrt(cov[1, 1]) * escala if con_color else 0.0
    return p[0], p[1], np.sqrt(cov[0, 0]) * escala, ec, rms, ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("carpeta", help="carpeta de resultados de la noche (con stack_B.fits y stack_V.fits)")
    ap.add_argument("--catalogo", choices=["apass", "gaia"], default="apass")
    ap.add_argument("--vmin", type=float, default=11.0, help="V mínima del catálogo (evita saturadas)")
    ap.add_argument("--vmax", type=float, default=16.0, help="V máxima del catálogo")
    ap.add_argument("--salida", default="resultados", help="carpeta donde guardar la tabla")
    ap.add_argument("--ra", help="AR del objeto (hh:mm:ss); si no, se busca por el nombre de OBJECT")
    ap.add_argument("--dec", help="Dec del objeto (dd:mm:ss)")
    ap.add_argument("--sin-color", action="store_true",
                    help="ajustar solo el ZP (término de color = 0); más robusto con pocas estrellas o poco rango de color")
    args = ap.parse_args()

    stacks = {f: leer_stack(args.carpeta, f) for f in ("B", "V")}
    dataV, satV, hdrV = stacks["V"]
    wcs0 = WCS(hdrV)
    ny, nx = dataV.shape
    # El centro de la imagen es el objeto (se apuntó a él). No se usa OBJCTRA/OBJCTDEC del header:
    # el 15/10/2026 quedaron ~40′ corridos respecto de NGC 7793.
    if args.ra and args.dec:
        obj = SkyCoord(args.ra, args.dec, unit=(u.hourangle, u.deg))
    else:
        obj = SkyCoord.from_name(str(hdrV.get("OBJECT", "")).strip())
    centro = obj
    print(f"Centro de búsqueda: {hdrV.get('OBJECT')} en {obj.to_string('hmsdms', precision=1)}")

    det = detectar(dataV, fwhm_px=12.0)
    det = det[:400]
    x, y = np.asarray(det["xcentroid"]), np.asarray(det["ycentroid"])
    print(f"Detectadas {len(det)} estrellas en V (se usan las 400 más brillantes)")

    radio = ESCALA * max(nx, ny) / 60 / np.sqrt(2) + 2
    cat = catalogo_apass(centro, radio) if args.catalogo == "apass" else catalogo_gaia(centro, radio)
    print(f"Catálogo {args.catalogo.upper()}: {len(cat)} estrellas con B y V en {radio:.0f}′")
    wcs, idx, ok = resolver_astrometria(wcs0, x, y, cat, centro)

    # Fotometría de apertura en ambos apilados (están en la misma grilla)
    fwhm = 12.0
    r_ap, r_in, r_out = 2.0 * fwhm, 3.0 * fwhm, 4.5 * fwhm
    m = {}
    for f in ("B", "V"):
        data, sat, hdr = stacks[f]
        flux, err = fotometria(data, x, y, r_ap, r_in, r_out)
        cerca_sat = np.array([sat[max(int(yy) - 10, 0):int(yy) + 11, max(int(xx) - 10, 0):int(xx) + 11].any()
                              for xx, yy in zip(x, y)])
        with np.errstate(invalid="ignore", divide="ignore"):
            m[f] = (-2.5 * np.log10(flux), 1.0857 * err / flux, cerca_sat, float(hdr.get("AIRMASS", np.nan)))

    pos = wcs.pixel_to_world(x, y)
    lejos = pos.separation(obj).arcmin > EXCLUIR_RADIO_ARCMIN
    borde = (x > r_out) & (x < nx - r_out) & (y > r_out) & (y < ny - r_out)
    c = cat[idx]
    sel = ok & lejos & borde & (c["V"] > args.vmin) & (c["V"] < args.vmax) & ~m["B"][2] & ~m["V"][2]
    color = np.asarray(c["B"] - c["V"])
    print(f"Estrellas útiles: {sel.sum()} (coinciden, fuera de {EXCLUIR_RADIO_ARCMIN:.0f}′ del objeto, "
          f"{args.vmin} < V < {args.vmax}, sin saturar)")

    filas = []
    for f in ("B", "V"):
        minst, einst, _, x_air = m[f]
        dm = np.where(sel, np.asarray(c[f]) - minst, np.nan)
        e = np.sqrt(einst ** 2 + np.nan_to_num(np.asarray(c["e" + f]), nan=0.05) ** 2)
        zp, cc, ezp, ecc, rms, usados = ajustar(dm, color, e, con_color=not args.sin_color)
        zp_fuera = zp + K[f] * x_air
        filas.append((f, zp, ezp, cc, ecc, rms, int(usados.sum()), x_air, zp_fuera, zp_fuera - ZP_2025[f]))
        print(f"\n{f}: ZP = {zp:.3f} ± {ezp:.3f} ADU/s a X = {x_air:.3f}  |  término de color c = {cc:+.3f} ± {ecc:.3f}"
              f"  |  rms {rms:.3f} mag con {usados.sum()} estrellas")
        print(f"   fuera de la atmósfera (k = {K[f]} supuesto): {zp_fuera:.3f}  (2025: {ZP_2025[f]:.3f}; "
              f"diferencia {zp_fuera - ZP_2025[f]:+.2f} mag)")

    salida = Path(args.salida)
    salida.mkdir(parents=True, exist_ok=True)
    noche = Path(args.carpeta).name
    res = Table(rows=filas, names=("filtro", "zp", "zp_err", "termino_color", "termino_color_err", "rms",
                                   "n_estrellas", "airmass", "zp_fuera_atmosfera", "dif_vs_2025"))
    sufijo = f"{args.catalogo}_sin_color" if args.sin_color else args.catalogo
    res.write(salida / f"zp_campo_{noche}_{sufijo}.csv", overwrite=True)
    est = Table({"x": x[sel], "y": y[sel], "ra": pos.ra.deg[sel], "dec": pos.dec.deg[sel],
                 "B_cat": c["B"][sel], "V_cat": c["V"][sel],
                 "b_inst": m["B"][0][sel], "v_inst": m["V"][0][sel]})
    est["B_V_cat"] = est["B_cat"] - est["V_cat"]
    est["dm_B"] = est["B_cat"] - est["b_inst"]
    est["dm_V"] = est["V_cat"] - est["v_inst"]
    est.write(salida / f"estrellas_campo_{noche}_{sufijo}.csv", overwrite=True)
    print("\nPor estrella (catálogo − instrumental):")
    for r in est:
        print(f"  V={r['V_cat']:6.2f}  B-V={r['B_V_cat']:5.2f}  dm_B={r['dm_B']:6.2f}  dm_V={r['dm_V']:6.2f}")
    print(f"\nTablas en {salida.resolve()}")


if __name__ == "__main__":
    sys.exit(main())
