"""Planificación: curvas de visibilidad, distancia a la Luna y bitácora de observación."""

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, get_body, get_sun
from astropy.table import Table
from astropy.time import Time
from astropy.utils import iers

from .io import exptime, filter_name, read_fits


def visibility(coord, location, date, utc_offset_h=-3.0, start_h=18.0, end_h=30.0, step_min=10):
    """Altura del objeto, del Sol y de la Luna durante la noche que empieza en ``date``.

    ``date`` es la fecha local (``"2026-10-08"``); las horas van de ``start_h``
    a ``end_h`` en hora local (30 = 06:00 del día siguiente).
    """
    hours = np.arange(start_h, end_h + 1e-9, step_min / 60.0)
    t0 = Time(f"{date}T00:00:00", scale="utc") - utc_offset_h * u.hour
    times = t0 + hours * u.hour
    with iers.conf.set_temp("auto_download", False):
        frame = AltAz(obstime=times, location=location)
        alt = coord.transform_to(frame)
        sun = get_sun(times).transform_to(frame)
        moon_gcrs = get_body("moon", times, location)
        moon = moon_gcrs.transform_to(frame)
        sep = moon_gcrs.separation(coord.transform_to(moon_gcrs.frame))
        sun_gcrs = get_body("sun", times, location)
        elong = sun_gcrs.separation(moon_gcrs)
    illum = (1 - np.cos(elong.rad)) / 2  # fracción iluminada (aprox. sin paralaje)
    tbl = Table()
    tbl["hora_local"] = hours
    tbl["alt_objeto"] = alt.alt.deg
    tbl["airmass"] = np.where(alt.alt.deg > 5, 1 / np.sin(np.radians(np.maximum(alt.alt.deg, 5))), np.nan)
    tbl["alt_sol"] = sun.alt.deg
    tbl["alt_luna"] = moon.alt.deg
    tbl["dist_luna"] = sep.deg
    tbl["iluminacion_luna"] = illum
    return tbl


def plot_visibility(ax, tbl, name, limit_hour=25.0):
    """Grafica la curva de visibilidad (horas > 24 son de la madrugada siguiente)."""
    h = np.asarray(tbl["hora_local"])
    ax.plot(h, tbl["alt_objeto"], "k-", label=name)
    ax.plot(h, tbl["alt_luna"], "--", color="0.5",
            label=f"Luna ({100 * np.nanmean(tbl['iluminacion_luna']):.0f} % iluminada)")
    dark = np.asarray(tbl["alt_sol"]) < -18
    ax.fill_between(h, 0, 90, where=dark, color="0.9", step="mid", label="noche astronómica")
    if limit_hour is not None:
        ax.axvline(limit_hour, color="r", lw=1, ls=":", label="límite 01:00")
    ax.axhline(30, color="0.6", lw=0.8)  # altura 30° ~ masa de aire 2
    ticks = np.arange(np.ceil(h[0]), h[-1] + 1, 2)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{int(t % 24):02d}:00" for t in ticks])
    ax.set_ylim(0, 90)
    ax.set_xlabel("Hora local")
    ax.set_ylabel("Altura [°]")
    ax.legend(fontsize=8, loc="upper right")
    return ax


def observation_log(files):
    """Tabla con los datos de cada archivo (para la sección "Observaciones" del informe)."""
    rows = []
    for f in files:
        _, h = read_fits(f)
        rows.append((f.name, str(h.get("IMAGETYP", "")), str(h.get("OBJECT", "")), filter_name(h) or "",
                     str(h.get("DATE-OBS", "")), exptime(h), float(h.get("AIRMASS", np.nan) or np.nan)))
    return Table(rows=rows, names=("archivo", "tipo", "objeto", "filtro", "date_obs_utc", "exptime", "airmass"))
