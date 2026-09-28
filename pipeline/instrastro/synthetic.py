"""Set de datos sintético para probar el pipeline sin tener observaciones.

Simula una cámara con bias con patrón, corriente oscura con píxeles
calientes, viñeteo y variaciones píxel a píxel del flat, ruido de Poisson y
de lectura, deriva del telescopio entre exposiciones, un cúmulo de estrellas
con colores conocidos sobre una nebulosa, y dos estrellas estándar. Como se
conocen los valores verdaderos (zeropoint, seeing, magnitudes), las pruebas
comprueban que el pipeline los recupera.
"""

from pathlib import Path

import numpy as np
import yaml
from astropy import units as u
from astropy.io import fits
from astropy.time import Time

from .photometry import tycho_to_johnson

TRUTH = {
    "zp": {"B": 20.8, "V": 21.2},
    "k": {"B": 0.25, "V": 0.15},
    "gain": 1.0,
    "readnoise": 5.0,
    "bias": 1000.0,
    "dark_current": 0.05,
    "saturation": 60000.0,
    "pixscale": 0.36,
}

STANDARDS = [
    {"nombre": "HIP 116375", "ra": "23:36:18", "dec": "-34:56:37", "BT": 8.453, "VT": 7.191},
    {"nombre": "HIP 117678", "ra": "23:53:17", "dec": "-32:58:44", "BT": 8.290, "VT": 7.160},
]

NIGHTS = [
    {"nombre": "noche1", "fecha": "2026-10-08", "fwhm_px": 7.0, "sky": {"B": 8.0, "V": 15.0},
     "utc": "2026-10-09T01:30:00", "airmass": 1.10},
    {"nombre": "noche2", "fecha": "2026-10-15", "fwhm_px": 9.0, "sky": {"B": 30.0, "V": 40.0},
     "utc": "2026-10-16T01:30:00", "airmass": 1.25},
]


def _render(img, x, y, flux, sigma):
    """Suma una gaussiana normalizada de flujo total ``flux`` centrada en (x, y)."""
    half = int(np.ceil(5 * sigma))
    xi, yi = int(round(x)), int(round(y))
    y0, y1 = max(yi - half, 0), min(yi + half + 1, img.shape[0])
    x0, x1 = max(xi - half, 0), min(xi + half + 1, img.shape[1])
    if y0 >= y1 or x0 >= x1:
        return
    yy, xx = np.mgrid[y0:y1, x0:x1]
    g = np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / sigma ** 2)
    img[y0:y1, x0:x1] += flux * g / (2 * np.pi * sigma ** 2)


def _write(path, data, **header):
    hdr = fits.Header()
    for key, value in header.items():
        hdr[key.replace("_", "-") if key == "DATE_OBS" else key] = value
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(np.clip(data, 0, 65535).astype(np.uint16), header=hdr).writeto(path, overwrite=True)


class _Camera:
    def __init__(self, size, rng):
        self.size = size
        self.rng = rng
        yy, xx = np.mgrid[:size, :size]
        r2 = ((xx - size / 2) ** 2 + (yy - size / 2) ** 2) / (size / 2) ** 2 / 2
        self.flat = {f: (1 - 0.25 * r2) * rng.normal(1, 0.01, (size, size)) for f in ("B", "V")}
        self.bias = TRUTH["bias"] + 3 * np.sin(xx / 40.0) + rng.normal(0, 1, size)[None, :]
        self.dark = np.full((size, size), TRUTH["dark_current"])
        hot = rng.integers(0, size, (2, 30))
        self.dark[hot[0], hot[1]] = 3.0

    def expose(self, rate, t, filt=None):
        """Imagen cruda en ADU a partir de un mapa de tasa (ADU/s) fuera de la atmósfera y del flat."""
        signal = self.dark * t
        if rate is not None:
            signal = signal + rate * self.flat[filt] * t
        electrons = self.rng.poisson(np.maximum(signal, 0) * TRUTH["gain"])
        return self.bias + electrons / TRUTH["gain"] + self.rng.normal(0, TRUTH["readnoise"], signal.shape)


def make_dataset(outdir, size=640, n_frames=5, seed=0, nights=NIGHTS):
    """Genera los FITS y un ``config.yaml`` en ``outdir``. Devuelve ``(verdad, ruta_config)``."""
    outdir = Path(outdir)
    rng = np.random.default_rng(seed)
    cam = _Camera(size, rng)
    zp, k = TRUTH["zp"], TRUTH["k"]

    # Cúmulo: estrellas con V y B-V conocidos, en coordenadas de la primera imagen V
    n_stars = 45
    margin = 40
    stars = {
        "x": rng.uniform(margin, size - margin, n_stars),
        "y": rng.uniform(margin, size - margin, n_stars),
        "V": rng.uniform(12.0, 16.5, n_stars),
        "BV": rng.uniform(0.1, 1.4, n_stars),
    }
    # Una estrella saturada para comprobar que se marca
    stars["x"][0], stars["y"][0], stars["V"][0], stars["BV"][0] = size * 0.7, size * 0.3, 7.5, 0.5
    stars["B"] = stars["V"] + stars["BV"]
    yy, xx = np.mgrid[:size, :size]
    nebula_amp = {"B": 10.0, "V": 20.0}

    truth = {**TRUTH, "stars": stars, "nights": {}, "size": size}
    night_cfgs = []
    bpm = np.ones((size, size), dtype=np.uint8)
    bpm[:, 100] = 0  # columna mala
    _write(outdir / "datos" / "BPM.fit", bpm)

    for n_i, night in enumerate(nights):
        d = outdir / "datos" / night["nombre"]
        sigma = night["fwhm_px"] / 2.3548
        t0 = Time(night["utc"])
        x_air = night["airmass"]
        common = {"OBJCTRA": "18 20 26", "OBJCTDEC": "-16 10 36"}

        for i in range(n_frames):
            _write(d / f"Bias_1x1_{i:04d}.fit", cam.expose(None, 0.0), IMAGETYP="Bias Frame", EXPTIME=0.0,
                   DATE_OBS=(t0 - 0.2 * u.day).isot)
            _write(d / f"Dark_100.000secs_{i:04d}.fit", cam.expose(None, 100.0), IMAGETYP="Dark Frame",
                   EXPTIME=100.0, DATE_OBS=(t0 - 0.2 * u.day).isot)
            for filt in ("B", "V"):
                t = 4.0 + i  # flats de crepúsculo con distinto tiempo y nivel
                level = rng.uniform(15000, 30000) / t
                _write(d / f"FLAT_{filt}_{i:03d}.fit", cam.expose(np.full((size, size), level), t, filt),
                       IMAGETYP="Flat Field", FILTER=filt, EXPTIME=t, DATE_OBS=(t0 - 0.1 * u.day).isot)

        shifts = {}
        for filt in ("V", "B"):
            for i in range(n_frames):
                dy, dx = (0.0, 0.0) if (filt == "V" and i == 0) else rng.uniform(-4, 4, 2)
                shifts[(filt, i)] = (dy, dx)
                rate = np.full((size, size), night["sky"][filt], dtype=float)
                rate += nebula_amp[filt] * np.exp(-((xx - size * 0.4 - dx) ** 2 + (yy - size * 0.55 - dy) ** 2)
                                                  / (2 * 60.0 ** 2))
                atten = 10 ** (-0.4 * k[filt] * x_air)
                for x, y, m in zip(stars["x"], stars["y"], stars[filt]):
                    _render(rate, x + dx, y + dy, 10 ** (-0.4 * (m - zp[filt])) * atten, sigma)
                t_obs = t0 + (i * 35 + (0 if filt == "V" else 200)) * u.s
                _write(d / f"OBJ_{filt}_30.000secs_{i:04d}.fit", cam.expose(rate, 30.0, filt), IMAGETYP="Light Frame",
                       OBJECT="NGC 6618", FILTER=filt, EXPTIME=30.0, DATE_OBS=t_obs.isot, AIRMASS=x_air, **common)
        truth["nights"][night["nombre"]] = {"fwhm_px": night["fwhm_px"], "shifts": shifts, "airmass": x_air}

        std_cfgs = []
        if n_i == 0:  # solo la primera noche tiene estándares (la segunda usa su ZP)
            for s_i, std in enumerate(STANDARDS):
                b_mag, v_mag = tycho_to_johnson(std["BT"], std["VT"])
                mags = {"B": b_mag, "V": v_mag}
                x_std = 1.05 + 0.3 * s_i
                files = {}
                for filt, t in (("B", 5.0), ("V", 3.0)):
                    for i in range(3):
                        rate = np.full((size, size), night["sky"][filt], dtype=float)
                        flux = 10 ** (-0.4 * (mags[filt] - zp[filt])) * 10 ** (-0.4 * k[filt] * x_std)
                        _render(rate, size / 2 + rng.uniform(-20, 20), size / 2 + rng.uniform(-20, 20), flux, sigma)
                        for _ in range(5):
                            _render(rate, *rng.uniform(40, size - 40, 2), 10 ** (-0.4 * (14 - zp[filt])), sigma)
                        name = f"HIP{s_i}_{filt}_{t:.3f}secs_{i:04d}.fit"
                        _write(d / "estandares" / name, cam.expose(rate, t, filt), IMAGETYP="Light Frame",
                               OBJECT=std["nombre"], FILTER=filt, EXPTIME=t, DATE_OBS=(t0 + 0.05 * u.day).isot,
                               AIRMASS=x_std, OBJCTRA=std["ra"].replace(":", " "),
                               OBJCTDEC=std["dec"].replace(":", " "))
                    files[filt] = f"datos/{night['nombre']}/estandares/HIP{s_i}_{filt}_*.fit"
                std_cfgs.append({"nombre": std["nombre"], "ra": std["ra"], "dec": std["dec"],
                                 "magnitudes": {"BT": std["BT"], "VT": std["VT"]}, "archivos": files})

        night_cfgs.append({
            "nombre": night["nombre"],
            "fecha": night["fecha"],
            "bias": f"datos/{night['nombre']}/Bias_*.fit",
            "darks": f"datos/{night['nombre']}/Dark_*.fit",
            "flats": {f: f"datos/{night['nombre']}/FLAT_{f}_*.fit" for f in ("B", "V")},
            "ciencia": {f: f"datos/{night['nombre']}/OBJ_{f}_*.fit" for f in ("B", "V")},
            "excluir": [],
            "wcs": {},
            "estandares": std_cfgs,
        })

    config = {
        "objeto": {"nombre": "NGC 6618 (sintético)", "ra": "18:20:26", "dec": "-16:10:36"},
        "sitio": {"nombre": "El Sauce", "latitud": -30.4725, "longitud": -70.7631, "altitud": 1600,
                  "utc_offset_h": -3},
        "instrumento": {"diametro_m": 0.5, "escala_arcsec_px": TRUTH["pixscale"], "ganancia_e_adu": TRUTH["gain"],
                        "saturacion_adu": TRUTH["saturation"], "bpm": "datos/BPM.fit",
                        "seeing_estimado_arcsec": 3.0},
        "filtros": ["V", "B"],
        "extincion": dict(k),
        "deteccion": {"extendida": {"box_size": 256, "n_pixels": 200}},
        "astrometria": {"modo": "header"},
        "salida": "resultados",
        "noches": night_cfgs,
    }
    config_path = outdir / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, allow_unicode=True, sort_keys=False), encoding="utf-8")
    return truth, config_path
