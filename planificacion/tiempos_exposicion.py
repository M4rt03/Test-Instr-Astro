"""Calculadora de tiempos de exposición para el MAS500 (Tarea 2, FIS2421).

Usa la ecuación de señal a ruido de un detector CCD/CMOS para estimar, en cada
filtro, cuánto debe durar cada exposición y cuántas tomar.

Solo usa la biblioteca estándar de Python:

    python tiempos_exposicion.py            # tablas en Markdown
    python tiempos_exposicion.py --zp-v 21.8 --ruido 3.2 --seeing 3.5

Los parámetros de INSTRUMENTO son estimaciones (ver tiempos_exposicion.md).
Conviene reemplazarlos por los medidos con los datos de 2025 (pipeline:
calibracion.json y zeropoints.ecsv) o con las exposiciones de prueba.
"""

import argparse
import math

# ----------------------------------------------------------------------------
# Instrumento y sitio
# ----------------------------------------------------------------------------
INSTRUMENTO = {
    "escala": 0.317,          # arcsec/px medida (medir_escala.py: 0.3168; astrometry.net 2025: 0.3163)
    "ganancia": 0.85,         # e-/ADU (header GAIN de los datos 2025)
    "ruido": 3.9,             # e- RMS, modo 16-bit HDR del GSENSE4040 FSI (ficha de Moravian)
    "oscuridad": 0.05,        # e-/s/px (supuesto: Moravian no lo publica; medir con los darks)
    "saturacion": 40000,      # ADU sobre el bias; el ADC llega a 65535 (pozo: 56600 e- = 66600 ADU)
    "lectura_s": 5,           # s entre exposiciones: descarga 0.25 s (USB 3) + guardado en TheSkyX (supuesto)
    # Zeropoints sobre la atmósfera: magnitud que da 1 e-/s. Estimados con
    # área efectiva 1570 cm2 (0.5 m con obstrucción ~45 %), flujo de Vega en B y V
    # (Bessell 1998) y eficiencia total 0.30 (B) y 0.48 (V).
    "zp_e": {"B": 21.98, "V": 22.05},
}
EXTINCION = {"B": 0.25, "V": 0.15}          # mag por masa de aire
# Brillo del cielo en el cenit (mag/arcsec2). Sin Luna: valores típicos de
# Cerro Tololo. Con Luna al 25 % (~5 días): valores típicos para esa fase.
CIELO = {
    "8/10 sin Luna": {"B": 22.7, "V": 21.8},
    "15/10 Luna 25 %": {"B": 22.1, "V": 21.6},
}
SEEING_SN = 4.0      # arcsec, para la S/N (informe 2025: 3.9-4.4")
SEEING_SAT = 2.5     # arcsec, para la saturación (caso con buen seeing = pico más alto)

# ----------------------------------------------------------------------------
# Objetos (candidatos_2026.md): V, B, tamaño (arcmin), altura de trabajo (°)
# ----------------------------------------------------------------------------
OBJETOS = [
    ("NGC 7293", 7.3, 7.5, (16.3, 16.3), 75),
    ("NGC 7009", 8.0, 8.3, (0.7, 0.7), 70),
    ("NGC 300", 8.7, 8.8, (19.4, 13.1), 70),
    ("NGC 7793", 9.3, 9.7, (10.4, 6.0), 75),
    ("NGC 247", 9.2, 9.7, (19.7, 5.5), 70),
    ("NGC 1097", 9.8, 10.1, (10.6, 6.4), 60),
    ("NGC 1291", 8.7, 9.4, (11.2, 9.9), 55),
    ("NGC 1316", 8.5, 9.4, (13.5, 7.7), 55),
    ("NGC 1313", 9.5, 9.7, (11.1, 9.1), 45),
    ("NGC 6744", 9.3, 9.1, (15.7, 9.8), 50),
]
ESTANDARES = [
    ("HD 195500", 7.32, 7.38), ("HD 202941", 7.07, 7.07), ("HD 210300", 6.44, 6.59),
    ("HD 215863", 7.69, 7.84), ("HD 220881", 7.45, 7.74), ("HD 562", 7.66, 7.80),
    ("HD 8130", 7.45, 7.50), ("HD 12206", 6.79, 6.81),
]

# Plan propuesto por objeto: (t_exp B, N B, t_exp V, N V)
PLAN = {
    "NGC 7009": (30, 12, 20, 12),     # núcleo muy brillante: exposiciones cortas
    "NGC 1097": (90, 9, 60, 9),       # núcleo Seyfert brillante
    "NGC 1316": (90, 9, 60, 9),       # centro muy concentrado
    "NGC 1291": (90, 9, 60, 9),
}
PLAN_DEFECTO = (120, 7, 90, 6)        # galaxias y nebulosas de bajo brillo superficial


def airmass(alt_deg):
    return 1.0 / math.sin(math.radians(alt_deg))


def zp_efectivo(filt, x, inst):
    """Magnitud que da 1 e-/s a masa de aire x."""
    return inst["zp_e"][filt] - EXTINCION[filt] * x


def tasa(m, zp):
    """e-/s de una fuente de magnitud m (o e-/s/arcsec2 si m es un brillo superficial)."""
    return 10 ** (0.4 * (zp - m))


def fraccion_pico(seeing, escala):
    """Fracción del flujo de una estrella gaussiana que cae en el píxel central."""
    fwhm_px = seeing / escala
    return 1.0 / (1.1331 * fwhm_px ** 2)


def snr(señal, cielo_px, npix, t, n, inst):
    """S/N del apilado de n exposiciones de t segundos (ecuación del CCD)."""
    s = señal * t
    ruido_px = cielo_px * t + inst["oscuridad"] * t + inst["ruido"] ** 2
    return n * s / math.sqrt(n * (s + npix * ruido_px))


def mu_limite(filt, x, mu_cielo, t, n, inst, objetivo=3.0):
    """Brillo superficial con S/N = objetivo en un elemento de seeing del apilado."""
    zp = zp_efectivo(filt, x, inst)
    area = math.pi / 4 * SEEING_SN ** 2
    npix = area / inst["escala"] ** 2
    cielo_px = tasa(mu_cielo, zp) * inst["escala"] ** 2
    lo, hi = 15.0, 32.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if snr(tasa(mid, zp) * area, cielo_px, npix, t, n, inst) > objetivo:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def tabla_objetos(inst, noche):
    filas = []
    for nombre, v, b, (a, c), alt in OBJETOS:
        x = airmass(alt)
        tb, nb, tv, nv = PLAN.get(nombre, PLAN_DEFECTO)
        area_obj = math.pi / 4 * a * c * 3600          # arcsec2
        fila = {"objeto": nombre, "X": x, "plan": (tb, nb, tv, nv)}
        for filt, m, t, n in (("B", b, tb, nb), ("V", v, tv, nv)):
            zp = zp_efectivo(filt, x, inst)
            mu = m + 2.5 * math.log10(area_obj)            # brillo superficial medio
            area_res = math.pi / 4 * SEEING_SN ** 2
            npix = area_res / inst["escala"] ** 2
            cielo_px = tasa(CIELO[noche][filt], zp) * inst["escala"] ** 2
            fila[filt] = {
                "mu": mu,
                "snr": snr(tasa(mu, zp) * area_res, cielo_px, npix, t, n, inst),
                "snr_borde": snr(tasa(mu + 2, zp) * area_res, cielo_px, npix, t, n, inst),
                "mu_lim": mu_limite(filt, x, CIELO[noche][filt], t, n, inst),
                "cielo_e": cielo_px * t,
                # brillo superficial que satura un píxel y estrella que satura
                "mu_sat": zp - 2.5 * math.log10(inst["saturacion"] * inst["ganancia"]
                                                / (t * inst["escala"] ** 2)),
                "m_sat": zp - 2.5 * math.log10(inst["saturacion"] * inst["ganancia"]
                                               / (t * fraccion_pico(SEEING_SAT, inst["escala"]))),
            }
        filas.append(fila)
    return filas


def tabla_estandares(inst, x=1.05):
    filas = []
    for nombre, v, b in ESTANDARES:
        fila = {"estrella": nombre}
        for filt, m in (("B", b), ("V", v)):
            zp = zp_efectivo(filt, x, inst)
            pico = tasa(m, zp) * fraccion_pico(SEEING_SAT, inst["escala"]) / inst["ganancia"]  # ADU/s
            t_sat = inst["saturacion"] / pico
            # tiempo redondeado hacia abajo a 0.5 s, con ~60 % del límite y mínimo 1 s
            t = max(1.0, math.floor(0.6 * t_sat * 2) / 2)
            señal = tasa(m, zp) * t
            npix = math.pi * (1.5 * SEEING_SN / inst["escala"]) ** 2   # apertura r = 1.5 FWHM
            ruido = math.sqrt(señal + npix * (inst["ruido"] ** 2 + tasa(CIELO["15/10 Luna 25 %"][filt], zp)
                                               * inst["escala"] ** 2 * t))
            fila[filt] = {"t_sat": t_sat, "t": t, "pico": pico * t, "snr": señal / ruido}
        filas.append(fila)
    return filas


def imprimir(inst):
    print("## Objetos: plan por bloque de 30 min\n")
    print("| Objeto | X | B: t_exp × N | V: t_exp × N | Tiempo total (min) | "
          "μB medio | μV medio | S/N B 8/10 · 15/10 | S/N V 8/10 · 15/10 | "
          "S/N B borde (μ+2) 8/10 | μ límite B 8/10 (S/N=3) | μ satura B · V | Estrella satura V |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    t8 = tabla_objetos(inst, "8/10 sin Luna")
    t15 = tabla_objetos(inst, "15/10 Luna 25 %")
    for f8, f15 in zip(t8, t15):
        tb, nb, tv, nv = f8["plan"]
        total = (nb * (tb + inst["lectura_s"]) + nv * (tv + inst["lectura_s"])) / 60
        print(f"| {f8['objeto']} | {f8['X']:.2f} | {tb} s × {nb} | {tv} s × {nv} | {total:.0f} | "
              f"{f8['B']['mu']:.1f} | {f8['V']['mu']:.1f} | "
              f"{f8['B']['snr']:.0f} · {f15['B']['snr']:.0f} | {f8['V']['snr']:.0f} · {f15['V']['snr']:.0f} | "
              f"{f8['B']['snr_borde']:.0f} | {f8['B']['mu_lim']:.1f} | "
              f"{f8['B']['mu_sat']:.1f} · {f8['V']['mu_sat']:.1f} | V < {f8['V']['m_sat']:.1f} |")

    print("\n## Estrellas estándar (X = 1.05, seeing 2.5\" para la saturación)\n")
    print("| Estrella | B | V | t satura B | t satura V | t_exp B | t_exp V | pico B (ADU) | pico V (ADU) | S/N B | S/N V |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for (nombre, v, b), f in zip(ESTANDARES, tabla_estandares(inst)):
        print(f"| {nombre} | {b:.2f} | {v:.2f} | {f['B']['t_sat']:.1f} s | {f['V']['t_sat']:.1f} s | "
              f"{f['B']['t']:.1f} s | {f['V']['t']:.1f} s | {f['B']['pico']:.0f} | {f['V']['pico']:.0f} | "
              f"{f['B']['snr']:.0f} | {f['V']['snr']:.0f} |")

    zp8 = {f: zp_efectivo(f, 1.2, inst) for f in "BV"}
    print("\n## Régimen de ruido (X = 1.2)\n")
    for noche, cielo in CIELO.items():
        for filt in "BV":
            cielo_px = tasa(cielo[filt], zp8[filt]) * inst["escala"] ** 2
            t_cielo = 10 * inst["ruido"] ** 2 / cielo_px
            print(f"- {noche}, {filt}: cielo {cielo_px:.3f} e-/s/px; el ruido del cielo domina "
                  f"(cielo > 10 × ruido²) recién con exposiciones de {t_cielo:.0f} s")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--zp-b", type=float, help="zeropoint B en e-/s sobre la atmósfera")
    p.add_argument("--zp-v", type=float, help="zeropoint V en e-/s sobre la atmósfera")
    p.add_argument("--ruido", type=float, help="ruido de lectura [e-]")
    p.add_argument("--oscuridad", type=float, help="corriente oscura [e-/s/px]")
    p.add_argument("--escala", type=float, help="escala de placa [arcsec/px]")
    p.add_argument("--seeing", type=float, help="seeing para la S/N [arcsec]")
    args = p.parse_args()

    global SEEING_SN
    inst = dict(INSTRUMENTO, zp_e=dict(INSTRUMENTO["zp_e"]))
    if args.zp_b is not None:
        inst["zp_e"]["B"] = args.zp_b
    if args.zp_v is not None:
        inst["zp_e"]["V"] = args.zp_v
    for clave in ("ruido", "oscuridad", "escala"):
        if getattr(args, clave) is not None:
            inst[clave] = getattr(args, clave)
    if args.seeing is not None:
        SEEING_SN = args.seeing
    imprimir(inst)


if __name__ == "__main__":
    main()
