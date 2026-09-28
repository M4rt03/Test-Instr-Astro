"""Orquestación del pipeline completo, noche por noche.

Cada función ``paso_*`` corresponde a un punto de la tarea y escribe sus
resultados (FITS, tablas ECSV/LaTeX, figuras PNG y un JSON de resumen) en
``<salida>/<noche>/``. ``run`` las ejecuta todas y al final compara las noches.
"""

import json
import logging
import warnings
from pathlib import Path

import numpy as np
from astropy.io import ascii
from astropy.stats import mad_std
from astropy.table import Table, vstack

from . import astrometry, color, detection, photometry, planning, plots, reduction, seeing
from .io import find_files, load_config, obs_time, read_fits, resolve, write_fits

log = logging.getLogger("instrastro")

# Radio (en FWHM) de la apertura "total" usada para las estándares y la corrección de apertura
APERTURE_FACTOR = 3.0


# ----------------------------------------------------------------------------
# Utilidades
# ----------------------------------------------------------------------------

class Context:
    """Configuración + rutas de una noche."""

    def __init__(self, cfg, night):
        self.cfg = cfg
        self.night = night
        self.name = night["nombre"]
        self.out = resolve(cfg, cfg.get("salida", "resultados")) / self.name
        self.figs = self.out / "figuras"
        self.figs.mkdir(parents=True, exist_ok=True)
        inst = cfg["instrumento"]
        self.pixscale = float(inst["escala_arcsec_px"])
        self.gain = float(inst.get("ganancia_e_adu", 1.0))
        self.saturation = float(inst.get("saturacion_adu", 60000))
        self.fwhm_guess = float(inst.get("seeing_estimado_arcsec", 3.0)) / self.pixscale
        self.filters = list(cfg.get("filtros", ["V", "B"]))
        self.location = photometry.site_location(cfg["sitio"])
        obj = cfg["objeto"]
        self.target = astrometry.target_coord(ra=obj.get("ra"), dec=obj.get("dec"))
        self.exclude = night.get("excluir", []) or []

    def files(self, pattern):
        return find_files(self.cfg, pattern, self.exclude)

    def k(self, filt):
        return float(self.cfg.get("extincion", {}).get(filt, 0.0))


def _save_json(path, data):
    def default(o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.bool_):
            return bool(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return str(o)
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False, default=default), encoding="utf-8")


def _save_table(tbl, path, latex=False):
    path = Path(path)
    tbl.write(path.with_suffix(".ecsv"), overwrite=True)
    if latex:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ascii.write(tbl, path.with_suffix(".tex"), format="latex", overwrite=True)


def _read_stack(path):
    from astropy.io import fits
    with fits.open(path) as hdul:
        data = hdul[0].data.astype(np.float32)
        header = hdul[0].header.copy()
        sat = hdul["SATMASK"].data.astype(bool) if "SATMASK" in hdul else np.zeros(data.shape, bool)
    return data, sat, header


# ----------------------------------------------------------------------------
# Pasos
# ----------------------------------------------------------------------------

def paso_bitacora(ctx):
    """Bitácora de todos los archivos de la noche y curvas de visibilidad."""
    night = ctx.night
    patterns = [night.get("bias"), night.get("darks")]
    patterns += list((night.get("flats") or {}).values()) + list((night.get("ciencia") or {}).values())
    for std in night.get("estandares") or []:
        patterns += list(std["archivos"].values())
    files = []
    for pat in patterns:
        files += ctx.files(pat)
    if files:
        _save_table(planning.observation_log(files), ctx.out / "bitacora", latex=True)

    if night.get("fecha"):
        targets = [(ctx.cfg["objeto"]["nombre"], ctx.target)]
        for std in night.get("estandares") or []:
            targets.append((std["nombre"], astrometry.target_coord(ra=std["ra"], dec=std["dec"])))
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, len(targets), figsize=(6 * len(targets), 4.5), squeeze=False)
        rows = []
        for ax, (name, coord) in zip(axes[0], targets):
            vis = planning.visibility(coord, ctx.location, night["fecha"], ctx.cfg["sitio"].get("utc_offset_h", -3))
            planning.plot_visibility(ax, vis, name)
            ax.set_title(f"{name} — noche del {night['fecha']}")
            best = int(np.nanargmax(vis["alt_objeto"]))
            rows.append((name, float(vis["hora_local"][best] % 24), float(vis["alt_objeto"][best]),
                         float(np.nanmin(vis["dist_luna"])), float(np.nanmean(vis["iluminacion_luna"]))))
        plots.save(fig, ctx.figs / "visibilidad.png")
        _save_table(Table(rows=rows, names=("objeto", "hora_culminacion", "alt_max", "dist_luna_min", "ilum_luna")),
                    ctx.out / "visibilidad", latex=True)


def paso_calibracion(ctx):
    """Master bias, dark current y master flats (punto 11: reducción)."""
    night = ctx.night
    bias_files = ctx.files(night["bias"])
    dark_files = ctx.files(night["darks"])
    if not bias_files or not dark_files:
        raise FileNotFoundError(f"[{ctx.name}] faltan bias o darks")
    log.info("[%s] master bias (%d) y dark current (%d)", ctx.name, len(bias_files), len(dark_files))
    bias = reduction.make_master_bias(bias_files)
    dark_current = reduction.make_dark_current(dark_files, bias)

    flats, levels = {}, {}
    for filt, pat in night["flats"].items():
        files = ctx.files(pat)
        log.info("[%s] master flat %s (%d)", ctx.name, filt, len(files))
        flats[filt], levels[filt] = reduction.make_master_flat(files, bias, dark_current, ctx.saturation)

    bad = None
    bpm_path = ctx.cfg["instrumento"].get("bpm")
    if bpm_path:
        bad = read_fits(resolve(ctx.cfg, bpm_path))[0] == 0

    # Ruido de lectura a partir de la diferencia de dos bias (elimina el patrón fijo)
    b1, b2 = read_fits(bias_files[0])[0], read_fits(bias_files[min(1, len(bias_files) - 1)])[0]
    readnoise = float(mad_std(b1 - b2) / np.sqrt(2)) if len(bias_files) > 1 else np.nan

    write_fits(ctx.out / "master_bias.fits", bias)
    write_fits(ctx.out / "dark_current.fits", dark_current)
    for filt, flat in flats.items():
        write_fits(ctx.out / f"master_flat_{filt}.fits", flat)
    plots.grid([bias, dark_current] + list(flats.values()),
               ["Master bias [ADU]", "Dark current [ADU/s]"] + [f"Master flat {f} (normalizado)" for f in flats],
               ctx.figs / "calibracion.png", label="")

    resumen = {
        "n_bias": len(bias_files), "n_darks": len(dark_files),
        "bias_mediana_adu": float(np.nanmedian(bias)),
        "ruido_lectura_adu": readnoise,
        "ruido_lectura_e": readnoise * ctx.gain,
        "dark_current_mediana_adu_s": float(np.nanmedian(dark_current)),
        "flats": {f: {"n": len(v), "niveles_adu": v} for f, v in levels.items()},
    }
    _save_json(ctx.out / "calibracion.json", resumen)
    return reduction.Calibration(bias=bias, dark_current=dark_current, flats=flats, bad_pixels=bad,
                                 saturation=ctx.saturation, flat_levels=levels)


def paso_ciencia(ctx, calib):
    """Reduce, alinea y apila las imágenes de ciencia de cada filtro en una grilla común."""
    def measure(red, sat):
        psf = seeing.measure_psf(red, sat, ctx.fwhm_guess, max_stars=60)
        s = seeing.summarize(psf, ctx.pixscale)
        return {"seeing_arcsec": s["seeing_arcsec"], "n_estrellas_psf": s["n_estrellas"],
                "elipticidad": s["elipticidad"]}

    def airmass(h):
        return photometry.airmass(h, ctx.location, ctx.target)

    stacks, reference, infos = {}, None, []
    for filt in ctx.filters:
        files = ctx.files(ctx.night["ciencia"][filt])
        if not files:
            raise FileNotFoundError(f"[{ctx.name}] no hay imágenes de ciencia en {filt}")
        log.info("[%s] reduciendo y apilando %d imágenes en %s", ctx.name, len(files), filt)
        stack, sat, header, info = reduction.reduce_and_stack(files, calib, filt, reference=reference,
                                                              measure=measure, airmass_func=airmass)
        if reference is None:
            reference = stack
        info["filtro"] = filt
        infos.append(info)
        write_fits(ctx.out / f"stack_{filt}.fits", stack, header, masks={"SATMASK": sat})
        stacks[filt] = (stack, sat, header)

        raw, raw_hdr = read_fits(files[0])
        red, _ = reduction.reduce_frame(raw, raw_hdr, calib, filt)
        plots.grid([raw, red, stack], [f"{files[0].name} (cruda) [ADU]", "reducida [ADU/s]",
                                       f"apilado de {len(files)} [ADU/s]"],
                   ctx.figs / f"reduccion_{filt}.png", label="")
    info = vstack(infos)
    _save_table(info, ctx.out / "exposiciones", latex=True)
    return stacks


def paso_astrometria(ctx, stacks):
    """Punto 12: solución WCS común a los apilados (están en la misma grilla)."""
    cfg_ast = ctx.cfg.get("astrometria", {}) or {}
    mode = cfg_ast.get("modo", "header")
    ref_filt = ctx.filters[0]
    data, _, header = stacks[ref_filt]

    if mode == "archivo":
        path = (ctx.night.get("wcs") or {}).get(ref_filt)
        if not path:
            raise ValueError(f"[{ctx.name}] modo 'archivo' requiere noches[].wcs.{ref_filt}")
        wcs = astrometry.load_wcs_file(resolve(ctx.cfg, path))
    elif mode == "astrometry_net":
        data_sub, _, rms = detection.subtract_background(data)
        params = {"n_sigma": 10.0, "n_pixels": 8, "kernel_fwhm_factor": 1.0, "deblend": False}
        segm, conv = detection.detect(data_sub, rms, ctx.fwhm_guess, params)
        src = detection.source_properties(data_sub, segm, conv)
        wcs = astrometry.solve_astrometry_net(src["x"], src["y"], src["flux"], data.shape, cfg_ast["api_key"],
                                              coord=ctx.target, pixscale_arcsec=ctx.pixscale)
    else:
        coord = ctx.target
        if "OBJCTRA" in header and "OBJCTDEC" in header:
            coord = astrometry.target_coord(header)
        wcs = astrometry.initial_wcs(data.shape, coord, ctx.pixscale)

    for filt, (stack, sat, hdr) in stacks.items():
        astrometry.apply_wcs(hdr, wcs)
        write_fits(ctx.out / f"stack_{filt}.fits", stack, hdr, masks={"SATMASK": sat})
    resumen = astrometry.wcs_summary(wcs, data.shape)
    resumen["modo"] = mode
    _save_json(ctx.out / "astrometria.json", resumen)
    return wcs


def paso_estandares(ctx, calib):
    """Punto 13: fotometría de apertura de las estrellas estándar, exposición por exposición."""
    rows = []
    for std in ctx.night.get("estandares") or []:
        mags = photometry.standard_magnitudes(std)
        coord = astrometry.target_coord(ra=std.get("ra"), dec=std.get("dec")) if std.get("ra") else None
        for filt, pat in std["archivos"].items():
            for f in ctx.files(pat):
                raw, hdr = read_fits(f)
                red, sat = reduction.reduce_frame(raw, hdr, calib, filt)
                m = photometry.measure_standard(red, sat, ctx.fwhm_guess, APERTURE_FACTOR)
                x_air = photometry.airmass(hdr, ctx.location, coord or astrometry.target_coord(hdr))
                zp = photometry.zeropoint(m["flux"], mags[filt], ctx.k(filt), x_air) if m["flux"] > 0 else np.nan
                rows.append((std["nombre"], filt, f.name, x_air, m["flux"], m["fwhm_px"] * ctx.pixscale,
                             m["r_ap"], m["saturated"], mags[filt], zp))
    if not rows:
        return None, None
    tbl = Table(rows=rows, names=("estrella", "filtro", "archivo", "airmass", "flujo_adu_s", "fwhm_arcsec",
                                  "r_apertura_px", "saturada", "m_std", "zp"))
    if np.any(tbl["saturada"]):
        log.warning("[%s] %d mediciones de estándares saturadas: se excluyen del ZP",
                    ctx.name, int(np.sum(tbl["saturada"])))
    resumen = photometry.summarize_zeropoints(tbl)
    _save_table(tbl, ctx.out / "estandares", latex=True)
    _save_table(resumen, ctx.out / "zeropoints", latex=True)
    return tbl, resumen


def paso_seeing(ctx, stacks):
    """Punto 15: seeing de cada exposición y variación espacial de la PSF en el apilado."""
    import matplotlib.pyplot as plt
    info = Table.read(ctx.out / "exposiciones.ecsv")
    D = float(ctx.cfg["instrumento"]["diametro_m"])
    resumen, psfs = {}, {}
    fig, axes = plt.subplots(2, len(stacks), figsize=(5.5 * len(stacks), 9), squeeze=False)
    for col, (filt, (stack, sat, _)) in enumerate(stacks.items()):
        psf = seeing.measure_psf(stack, sat, ctx.fwhm_guess)
        psfs[filt] = psf
        s = seeing.summarize(psf, ctx.pixscale)
        per_exp = np.asarray(info["seeing_arcsec"][info["filtro"] == filt], dtype=float)
        grid, counts = seeing.spatial_grid(psf, stack.shape)
        slope, slope_err = seeing.radial_trend(psf, stack.shape)
        diff = seeing.diffraction_limit(filt, D)
        resumen[filt] = {
            "seeing_apilado_arcsec": s["seeing_arcsec"],
            "dispersion_arcsec": s["dispersion_arcsec"],
            "elipticidad_mediana": s["elipticidad"],
            "n_estrellas": s["n_estrellas"],
            "seeing_por_exposicion_arcsec": per_exp.tolist(),
            "seeing_exposiciones_mediana": float(np.nanmedian(per_exp)) if per_exp.size else np.nan,
            "limite_difraccion": diff,
            "razon_seeing_difraccion": s["seeing_arcsec"] / diff["fwhm_airy_arcsec"],
            "fwhm_grilla_arcsec": (grid * ctx.pixscale).tolist(),
            "estrellas_por_celda": counts.tolist(),
            "pendiente_radial_arcsec_por_1000px": slope * ctx.pixscale,
            "pendiente_radial_err": slope_err * ctx.pixscale,
        }
        _save_table(psf, ctx.out / f"psf_{filt}")

        ax = axes[0, col]
        ax.hist(np.asarray(psf["fwhm_px"]) * ctx.pixscale, bins=25, color="0.6")
        ax.axvline(s["seeing_arcsec"], color="k", label=f"mediana {s['seeing_arcsec']:.2f}″")
        ax.axvline(diff["fwhm_airy_arcsec"], color="r", ls="--", label=f"difracción {diff['fwhm_airy_arcsec']:.2f}″")
        ax.set_xlabel("FWHM [arcsec]")
        ax.set_title(f"Filtro {filt}: {s['n_estrellas']} estrellas")
        ax.legend()
        ax = axes[1, col]
        im = ax.imshow(grid * ctx.pixscale, origin="lower", cmap="viridis",
                       extent=(0, stack.shape[1], 0, stack.shape[0]))
        for j in range(grid.shape[0]):
            for i in range(grid.shape[1]):
                if np.isfinite(grid[j, i]):
                    ax.text((i + 0.5) * stack.shape[1] / grid.shape[1], (j + 0.5) * stack.shape[0] / grid.shape[0],
                            f"{grid[j, i] * ctx.pixscale:.2f}″\n(n={counts[j, i]})", ha="center", va="center",
                            color="w", fontsize=9)
        ax.set_title(f"FWHM por zona del detector ({filt})")
        fig.colorbar(im, ax=ax, label="FWHM [arcsec]", fraction=0.046, pad=0.04)
    plots.save(fig, ctx.figs / "seeing.png")
    _save_json(ctx.out / "seeing.json", resumen)
    return resumen, psfs


def paso_catalogos(ctx, stacks, wcs, zps, psfs):
    """Punto 14: catálogos con dos configuraciones de detección y diagrama color-magnitud."""
    import matplotlib.pyplot as plt
    configs = {name: {**params, **((ctx.cfg.get("deteccion") or {}).get(name) or {})}
               for name, params in detection.DEFAULT_CONFIGS.items()}
    seeing_px = {f: float(np.median(psfs[f]["fwhm_px"])) if len(psfs[f]) else ctx.fwhm_guess for f in stacks}
    catalogs, resumen = {}, {}
    fig, axes = plt.subplots(len(configs), len(stacks), figsize=(6 * len(stacks), 6 * len(configs)), squeeze=False)
    for row, (name, params) in enumerate(configs.items()):
        for col, (filt, (stack, sat, header)) in enumerate(stacks.items()):
            log.info("[%s] catálogo %s en %s", ctx.name, name, filt)
            gain_eff = ctx.gain * float(header.get("TOTEXP", 1.0))
            tbl, segm, data_sub, background, rms = detection.build_catalog(
                stack, seeing_px[filt], params, gain_eff=gain_eff, sat_mask=sat, wcs=wcs)
            info = {"n_fuentes": len(tbl), "parametros": params}
            if segm is None:
                resumen[f"{name}_{filt}"] = info
                continue
            if params.get("flux") == "aperture":
                fwhm = seeing_px[filt]
                r_ap = params.get("r_apertura_fwhm", 1.0) * fwhm
                r_in, r_out = (APERTURE_FACTOR + 0.5) * fwhm, (APERTURE_FACTOR + 2.0) * fwhm
                # Se mide sobre la imagen sin fondo ni nebulosa (data_sub): el anillo local por sí solo
                # no sigue la curvatura de la nebulosa y sesga a las estrellas débiles que están sobre ella.
                flux, err, _ = photometry.aperture_photometry_list(data_sub, tbl["x"], tbl["y"], r_ap, r_in, r_out,
                                                                   gain_eff, float(header.get("NOISECOR", 1.0)))
                # Estrellas no saturadas del propio catálogo para construir la PSF empírica
                stars = ~np.asarray(tbl["saturated"]) & (np.asarray(tbl["elongation"]) < 1.5)
                corr, corr_err, n = photometry.aperture_correction(
                    data_sub, np.asarray(tbl["x"])[stars], np.asarray(tbl["y"])[stars], r_ap,
                    APERTURE_FACTOR * fwhm, r_in, r_out)
                if n < 5:
                    log.warning("[%s] solo %d estrellas para la corrección de apertura en %s", ctx.name, n, filt)
                # El error de la corrección es sistemático (igual para todas): se suma en cuadratura
                rel_err = corr_err / corr if np.isfinite(corr_err) else 0.0
                tbl["flux"] = flux * corr
                tbl["flux_err"] = np.hypot(err / flux, rel_err) * np.abs(tbl["flux"])
                info["radio_apertura_px"] = r_ap
                info["correccion_apertura_mag"] = float(-2.5 * np.log10(corr))
                info["correccion_apertura_err_mag"] = float(2.5 / np.log(10) * rel_err)
                info["n_estrellas_correccion"] = n
            zp = zps.get(filt)
            if zp is not None:
                detection.calibrate(tbl, zp["zp"], zp["zp_err"] if np.isfinite(zp["zp_err"]) else 0.0,
                                    ctx.k(filt), float(header.get("AIRMASS", 1.0)))
            info["area_mediana_px"] = float(np.median(tbl["area"])) if len(tbl) else np.nan
            info["area_max_px"] = float(np.max(tbl["area"])) if len(tbl) else np.nan
            catalogs[(name, filt)] = tbl
            resumen[f"{name}_{filt}"] = info
            _save_table(tbl, ctx.out / f"catalogo_{name}_{filt}")
            plots.segmentation(axes[row, col], segm, f"{name}, filtro {filt}")
    plots.save(fig, ctx.figs / "segmentacion.png")

    # Diagrama color-magnitud con el catálogo agresivo (el que separa estrellas)
    if ("agresiva", "B") in catalogs and ("agresiva", "V") in catalogs and "mag" in catalogs[("agresiva", "V")].colnames:
        tol = max(2.0, 0.7 * np.mean(list(seeing_px.values())))
        cmd = detection.crossmatch(catalogs[("agresiva", "B")], catalogs[("agresiva", "V")], tol)
        _save_table(cmd, ctx.out / "catalogo_BV")
        good = ~cmd["saturated"] & np.isfinite(cmd["B_V"]) & (cmd["B_V_err"] < 0.2)
        fig, ax = plt.subplots(figsize=(6, 7))
        ax.errorbar(cmd["B_V"][good], cmd["V"][good], xerr=cmd["B_V_err"][good], yerr=cmd["V_err"][good],
                    fmt="o", ms=2, color="k", ecolor="0.7", elinewidth=0.5)
        ax.invert_yaxis()
        ax.set_xlabel("B − V")
        ax.set_ylabel("V")
        ax.set_title(f"Diagrama color-magnitud ({int(good.sum())} estrellas, σ(B−V) < 0.2)")
        plots.save(fig, ctx.figs / "color_magnitud.png")
        resumen["n_cruzadas_BV"] = len(cmd)
        resumen["n_diagrama_cm"] = int(good.sum())
    _save_json(ctx.out / "catalogos.json", resumen)
    return catalogs, resumen


def paso_color(ctx, stacks, zps):
    """Punto 16: imagen RGB (R=V, G=(B+V)/2, B=B) y mapa de color B-V."""
    import matplotlib.pyplot as plt
    if "B" not in stacks or "V" not in stacks:
        return
    img_b, img_v = stacks["B"][0], stacks["V"][0]
    rgb = color.rgb_image(img_b, img_v)
    plt.imsave(ctx.figs / "color_rgb.png", rgb, origin="lower")
    if zps.get("B") is not None and zps.get("V") is not None:
        zp_b = zps["B"]["zp"] - ctx.k("B") * float(stacks["B"][2].get("AIRMASS", 1.0))
        zp_v = zps["V"]["zp"] - ctx.k("V") * float(stacks["V"][2].get("AIRMASS", 1.0))
        bv = color.color_map(img_b, img_v, zp_b, zp_v)
        fig, ax = plt.subplots(figsize=(7, 6))
        lo, hi = np.nanpercentile(bv, [5, 95]) if np.isfinite(bv).any() else (0, 1)
        im = ax.imshow(bv, origin="lower", cmap="RdYlBu_r", vmin=lo, vmax=hi)
        fig.colorbar(im, ax=ax, label="B − V [mag]")
        ax.set_title("Mapa de color B − V (suavizado)")
        plots.save(fig, ctx.figs / "mapa_color_BV.png")


def paso_resumen_noche(ctx, stacks, zps, seeing_res):
    """Indicadores para comparar noches: fondo de cielo, magnitud límite, Luna, seeing."""
    from astropy.coordinates import get_body
    from astropy.utils import iers

    rows = []
    for filt, (stack, _, header) in stacks.items():
        _, background, rms = detection.subtract_background(stack, box_size=256)
        sky = float(np.nanpercentile(background, 20))  # evita la nebulosa
        noise = float(np.nanmedian(rms)) * float(header.get("NOISECOR", 1.0))
        x_air = float(header.get("AIRMASS", np.nan))
        zp = zps.get(filt)
        fwhm_px = seeing_res[filt]["seeing_apilado_arcsec"] / ctx.pixscale
        sky_mag = mlim = np.nan
        if zp is not None and np.isfinite(zp["zp"]):
            zp_eff = zp["zp"] - ctx.k(filt) * x_air
            sky_mag = zp_eff - 2.5 * np.log10(sky / ctx.pixscale ** 2) if sky > 0 else np.nan
            mlim = detection.limiting_magnitude(noise, fwhm_px, zp_eff)
        t = obs_time(header)
        with iers.conf.set_temp("auto_download", False):
            moon = get_body("moon", t, ctx.location)
            sun = get_body("sun", t, ctx.location)
            sep = moon.separation(ctx.target.transform_to(moon.frame)).deg
            illum = (1 - np.cos(sun.separation(moon).rad)) / 2
        rows.append((ctx.name, filt, int(header.get("NCOMBINE", 1)), float(header.get("TOTEXP", np.nan)), x_air,
                     seeing_res[filt]["seeing_apilado_arcsec"], seeing_res[filt]["seeing_exposiciones_mediana"],
                     sky, sky_mag, noise, mlim, zp["zp"] if zp is not None else np.nan, float(illum), float(sep)))
    return Table(rows=rows, names=("noche", "filtro", "n_exp", "t_total_s", "airmass", "seeing_apilado",
                                   "seeing_exposiciones", "cielo_adu_s_px", "cielo_mag_arcsec2", "rms_adu_s",
                                   "mag_limite_5sigma", "zp", "luna_iluminacion", "luna_distancia_deg"))


# ----------------------------------------------------------------------------
# Ejecución completa
# ----------------------------------------------------------------------------

def _zp_dict(resumen):
    if resumen is None:
        return {}
    return {row["filtro"]: {"zp": float(row["zp"]), "zp_err": float(row["zp_err"])} for row in resumen
            if np.isfinite(row["zp"])}


def run(config_path, noches=None):
    """Ejecuta todo el pipeline. ``noches`` permite correr solo algunas (por nombre)."""
    cfg = load_config(config_path)
    nights = [n for n in cfg["noches"] if noches is None or n["nombre"] in noches]
    contexts, calibs, zps = [], {}, {}

    # 1) Calibraciones y estándares de todas las noches (el ZP puede compartirse)
    for night in nights:
        ctx = Context(cfg, night)
        contexts.append(ctx)
        paso_bitacora(ctx)
        calibs[ctx.name] = paso_calibracion(ctx)
        _, resumen = paso_estandares(ctx, calibs[ctx.name])
        zps[ctx.name] = _zp_dict(resumen)

    all_zp = [z for z in zps.values() if z]
    global_zp = {}
    for filt in {f for z in all_zp for f in z}:
        vals = [z[filt]["zp"] for z in all_zp if filt in z]
        errs = [z[filt]["zp_err"] for z in all_zp if filt in z]
        global_zp[filt] = {"zp": float(np.mean(vals)),
                           "zp_err": float(np.nanmax(errs + [np.std(vals) if len(vals) > 1 else 0.0]))}

    # 2) Ciencia, astrometría, seeing, catálogos y color por noche
    comparison = []
    for ctx in contexts:
        night_zp = zps[ctx.name] or global_zp
        if not zps[ctx.name]:
            log.warning("[%s] sin estándares propias: se usa el ZP promedio de las otras noches", ctx.name)
        stacks = paso_ciencia(ctx, calibs[ctx.name])
        wcs = paso_astrometria(ctx, stacks)
        seeing_res, psfs = paso_seeing(ctx, stacks)
        paso_catalogos(ctx, stacks, wcs, night_zp, psfs)
        paso_color(ctx, stacks, night_zp)
        comparison.append(paso_resumen_noche(ctx, stacks, night_zp, seeing_res))
        del stacks

    out = resolve(cfg, cfg.get("salida", "resultados"))
    if comparison:
        tbl = vstack(comparison)
        for col in tbl.colnames:
            if tbl[col].dtype.kind == "f":
                tbl[col].format = ".3f"
        _save_table(tbl, out / "comparacion_noches", latex=True)
        log.info("Comparación entre noches:\n%s", tbl)
    return out
