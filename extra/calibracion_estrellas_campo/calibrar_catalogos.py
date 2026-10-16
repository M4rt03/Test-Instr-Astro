"""Calibra los catálogos del pipeline con el ZP de estrellas del campo y dibuja el diagrama
color-magnitud (fuera de la tarea).

- V: ZP con las estrellas de APASS, sin término de color (no queda determinado en V).
- B: ZP con las estrellas de APASS y término de color (c ≈ +2, posible fuga al rojo del filtro B).
El ajuste se hace sobre las magnitudes instrumentales del propio catálogo (misma apertura).

Con esas ecuaciones,
    V       = v + ZP_V
    B − V   = [(b − v) + ZP_B − ZP_V] / (1 − c_B)
    B       = V + (B − V)
donde b y v son magnitudes instrumentales totales (catálogo "agresiva" del pipeline, corregido
por apertura). El color solo está calibrado en el rango de colores de las estrellas de la
calibración (B−V ≈ 0,45–1,12); fuera de ese rango es una extrapolación.

    python calibrar_catalogos.py ../../pipeline/resultados/noche_15oct
"""

import json
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from scipy.spatial import cKDTree

from calibrar_campo import ajustar

warnings.simplefilter("ignore")

# NGC 7793 (SIMBAD): centro, diámetros 10,0′ × 6,0′, ángulo de posición 84°
GAL = SkyCoord("23h57m49.75s -32d35m27.7s")
GAL_A, GAL_B, GAL_PA = 10.0 / 2, 6.0 / 2, 84.0
SN_MIN = 5.0
COLOR_CAL = (0.45, 1.12)      # rango de B−V de las estrellas usadas en la calibración


def main():
    carpeta = Path(sys.argv[1])
    aqui = Path(__file__).parent
    res = aqui / "resultados"
    noche = carpeta.name
    corr = json.load(open(carpeta / "catalogos.json"))

    cat = {}
    for f in ("B", "V"):
        t = Table.read(carpeta / f"catalogo_agresiva_{f}.ecsv")
        t = t[(t["flux"] > 0) & (t["flux"] / t["flux_err"] > SN_MIN) & ~t["saturated"]]
        ap = corr[f"agresiva_{f}"]["correccion_apertura_mag"] or 0.0
        t["m_inst"] = -2.5 * np.log10(t["flux"]) + ap
        t["e_inst"] = 1.0857 * t["flux_err"] / t["flux"]
        cat[f] = t

    # cruce B-V por posición: los apilados están en la misma grilla
    tb, tv = cat["B"], cat["V"]
    d, i = cKDTree(np.column_stack([tb["x"], tb["y"]])).query(np.column_stack([tv["x"], tv["y"]]))
    ok = d < 3.0
    tv, tb = tv[ok], tb[i[ok]]

    # El ZP se ajusta de nuevo sobre las magnitudes instrumentales de ESTE catálogo, con las mismas
    # estrellas de APASS: los catálogos del pipeline usan otra apertura (y su corrección) que la
    # fotometría de calibrar_campo.py, y mezclar ambos sistemas deja ~0,1 mag de diferencia.
    ref = Table.read(res / f"estrellas_campo_{noche}_apass.csv")
    dr, ir = cKDTree(np.column_stack([tv["x"], tv["y"]])).query(np.column_stack([ref["x"], ref["y"]]))
    m = dr < 5
    color_ref = np.asarray(ref["B_cat"] - ref["V_cat"])[m]
    ev = np.asarray(tv["e_inst"])[ir[m]]
    eb = np.asarray(tb["e_inst"])[ir[m]]
    zv, _, ezv, _, rms_v, n_v = ajustar(np.asarray(ref["V_cat"])[m] - np.asarray(tv["m_inst"])[ir[m]],
                                        color_ref, np.hypot(ev, 0.03), con_color=False)
    zb, cb, ezb, ecb, rms_b, n_b = ajustar(np.asarray(ref["B_cat"])[m] - np.asarray(tb["m_inst"])[ir[m]],
                                           color_ref, np.hypot(eb, 0.03))
    print(f"ZP sobre el catálogo del pipeline ({m.sum()} estrellas de APASS): "
          f"V = {zv:.3f} ± {ezv:.3f} (rms {rms_v:.3f}); "
          f"B = {zb:.3f} ± {ezb:.3f} con c_B = {cb:+.2f} ± {ecb:.2f} (rms {rms_b:.3f})")

    V = tv["m_inst"] + zv
    BV = ((tb["m_inst"] - tv["m_inst"]) + zb - zv) / (1 - cb)
    B = V + BV
    eV = np.hypot(tv["e_inst"], ezv)
    eBV = np.hypot(np.hypot(tb["e_inst"], tv["e_inst"]), np.hypot(ezb, ezv)) / abs(1 - cb)

    # ¿dentro del disco de NGC 7793? (elipse de 10′ × 6′)
    pos = SkyCoord(tv["ra"], tv["dec"], unit="deg")
    off = pos.transform_to(GAL.skyoffset_frame())
    xi, eta = off.lon.arcmin, off.lat.arcmin
    pa = np.radians(GAL_PA)
    mayor = xi * np.sin(pa) + eta * np.cos(pa)
    menor = xi * np.cos(pa) - eta * np.sin(pa)
    disco = (mayor / GAL_A) ** 2 + (menor / GAL_B) ** 2 <= 1

    out = Table({"ra": tv["ra"], "dec": tv["dec"], "x": tv["x"], "y": tv["y"],
                 "V": V, "e_V": eV, "B": B, "B_V": BV, "e_B_V": eBV,
                 "en_disco": disco, "elongacion": tv["elongation"],
                 "color_en_rango_calibrado": (BV >= COLOR_CAL[0]) & (BV <= COLOR_CAL[1])})
    out.write(res / f"catalogo_calibrado_{noche}.csv", overwrite=True)
    print(f"{len(out)} fuentes en B y V con S/N > {SN_MIN:.0f}: {disco.sum()} en el disco de NGC 7793, "
          f"{(~disco).sum()} fuera")
    print(f"V: {np.min(V):.2f} a {np.max(V):.2f}   B−V: mediana {np.median(BV):.2f} "
          f"(disco {np.median(BV[disco]):.2f}, campo {np.median(BV[~disco]):.2f})")
    print(f"En el rango de color calibrado ({COLOR_CAL[0]}–{COLOR_CAL[1]}): "
          f"{out['color_en_rango_calibrado'].sum()} de {len(out)}")

    # --- Diagrama color-magnitud
    tinta, tinta2, grilla = "#0b0b0b", "#52514e", "#e4e3df"
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": tinta2, "axes.labelcolor": tinta,
                         "xtick.color": tinta2, "ytick.color": tinta2})
    fig, ax = plt.subplots(figsize=(6.4, 6.0), dpi=150)
    fig.patch.set_facecolor("#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    ax.axvspan(*COLOR_CAL, color="#ecebe7", zorder=0, lw=0)
    ax.text(np.mean(COLOR_CAL), 0.985, "color calibrado", transform=ax.get_xaxis_transform(),
            ha="center", va="top", color=tinta2, fontsize=8)
    for sel, color, nombre in ((~disco, "#2a78d6", f"Campo (estrellas de la Vía Láctea), N = {(~disco).sum()}"),
                               (disco, "#eb6834", f"Disco de NGC 7793, N = {disco.sum()}")):
        ax.errorbar(BV[sel], V[sel], xerr=eBV[sel], fmt="none", ecolor=color, alpha=0.25, lw=0.8, zorder=2)
        ax.scatter(BV[sel], V[sel], s=18, color=color, edgecolor="#fcfcfb", linewidth=0.6, label=nombre, zorder=3)
    ax.invert_yaxis()
    ax.set_xlabel("B − V (mag)")
    ax.set_ylabel("V (mag)")
    ax.set_title("NGC 7793, 15/10/2026: diagrama color-magnitud\n"
                 "ZP con estrellas del campo (APASS); color con término de color c_B = "
                 f"{cb:+.2f}", fontsize=10, color=tinta, loc="left")
    ax.grid(color=grilla, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, loc="upper right", fontsize=8.5, labelcolor=tinta)
    fig.tight_layout()
    salida = res / f"cmd_ngc7793_{noche}.png"
    fig.savefig(salida, facecolor=fig.get_facecolor())
    print(f"Diagrama en {salida}")


if __name__ == "__main__":
    main()
