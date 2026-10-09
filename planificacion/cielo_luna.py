"""Brillo del cielo con Luna y S/N de los objetos del grupo en las noches del 15/10 y 22/10.

Calcula, hora por hora, la posición de la Luna (fórmulas de baja precisión del Astronomical
Almanac, ~0,3°), su altura y su fase, la distancia a cada objeto y el brillo del cielo en V con
el modelo de Krisciunas & Schaefer (1991, PASP 103, 1033). Para B, la luz de Luna dispersada se
supone con B−V = 0,0 ± 0,3 (luz solar enrojecida por la Luna y azulada por la dispersión de
Rayleigh); por eso μB es más incierto que μV (±0,3–0,5 mag).

Con ese cielo recalcula la S/N del plan de tiempos_exposicion.py.

Solo usa la biblioteca estándar:

    python cielo_luna.py
"""

import math

import tiempos_exposicion as te

D2R = math.pi / 180
SITIO = (-30.4597, -70.7503)          # El Sauce (header de 2025)
UTC_LOCAL = -3                         # hora de Chile en octubre
CIELO_OSCURO = {"B": 22.7, "V": 21.8}  # cenit sin Luna, mag/arcsec2 (igual que tiempos_exposicion.py)
COLOR_LUZ_LUNA = 0.0                   # B−V supuesto del cielo iluminado por la Luna

OBJETOS = {  # nombre: (AR, Dec en grados)
    "NGC 7793": (359.4573, -32.5910),
    "NGC 1291": (49.3274, -41.1080),
}
NOCHES = {"15/10": (2026, 10, 15), "22/10": (2026, 10, 22)}
HORAS = [21, 22, 23, 24, 25, 26]       # hora local (24 = 00:00 del día siguiente)


def jd_local(fecha, hora):
    a, m, d = fecha
    # JD a las 0 h UT de la fecha (algoritmo de Meeus para fechas gregorianas)
    y, mm = (a, m) if m > 2 else (a - 1, m + 12)
    jd0 = int(365.25 * (y + 4716)) + int(30.6001 * (mm + 1)) + d + 2 - y // 100 + y // 400 - 1524.5
    return jd0 + (hora - UTC_LOCAL) / 24


def ecl_a_ecuatorial(lam, beta, T):
    eps = (23.439 - 0.013 * T) * D2R
    lam, beta = lam * D2R, beta * D2R
    ra = math.atan2(math.sin(lam) * math.cos(eps) - math.tan(beta) * math.sin(eps), math.cos(lam))
    dec = math.asin(math.sin(beta) * math.cos(eps) + math.cos(beta) * math.sin(eps) * math.sin(lam))
    return ra / D2R % 360, dec / D2R


def sol(jd):
    n = jd - 2451545.0
    L = (280.460 + 0.9856474 * n) % 360
    g = (357.528 + 0.9856003 * n) * D2R
    return (L + 1.915 * math.sin(g) + 0.020 * math.sin(2 * g)) % 360     # longitud eclíptica


def luna(jd):
    """Longitud y latitud eclípticas y paralaje horizontal de la Luna, en grados."""
    T = (jd - 2451545.0) / 36525
    s = lambda a, b: math.sin((a + b * T) * D2R)
    c = lambda a, b: math.cos((a + b * T) * D2R)
    lam = (218.32 + 481267.881 * T + 6.29 * s(135.0, 477198.87) - 1.27 * s(259.3, -413335.36)
           + 0.66 * s(235.7, 890534.22) + 0.21 * s(269.9, 954397.74) - 0.19 * s(357.5, 35999.05)
           - 0.11 * s(186.5, 966404.03)) % 360
    beta = (5.13 * s(93.3, 483202.02) + 0.28 * s(228.2, 960400.89) - 0.28 * s(318.3, 6003.15)
            - 0.17 * s(217.6, -407332.21))
    paralaje = (0.9508 + 0.0518 * c(135.0, 477198.87) + 0.0095 * c(259.3, -413335.36)
                + 0.0078 * c(235.7, 890534.22) + 0.0028 * c(269.9, 954397.74))
    return lam, beta, paralaje, T


def altura(ra, dec, jd):
    lst = (280.46061837 + 360.98564736629 * (jd - 2451545.0) + SITIO[1]) % 360
    h = (lst - ra) * D2R
    lat, dec = SITIO[0] * D2R, dec * D2R
    return math.asin(math.sin(lat) * math.sin(dec) + math.cos(lat) * math.cos(dec) * math.cos(h)) / D2R


def separacion(ra1, de1, ra2, de2):
    c = (math.sin(de1 * D2R) * math.sin(de2 * D2R)
         + math.cos(de1 * D2R) * math.cos(de2 * D2R) * math.cos((ra1 - ra2) * D2R))
    return math.acos(max(-1.0, min(1.0, c))) / D2R


def x_ks(z_deg):
    """Masa de aire del modelo de Krisciunas & Schaefer."""
    return (1 - 0.96 * math.sin(z_deg * D2R) ** 2) ** -0.5


def nl_a_mag(b_nl):
    return (20.7233 - math.log(b_nl / 34.08)) / 0.92104



def brillo_cielo(alt_obj, alt_luna, rho, fase_deg):
    """μB y μV (mag/arcsec2) en la dirección del objeto.

    Cielo oscuro: el del cenit, más brillante hacia el horizonte según Krisciunas & Schaefer.
    Luz de Luna: modelo de K&S en V; en B se usa la misma luz con el color COLOR_LUZ_LUNA.
    """
    z, zm = 90 - alt_obj, 90 - alt_luna
    luz_v = 0.0                                                  # nanoLamberts, en V
    if alt_luna > 0:
        kv = te.EXTINCION["V"]
        i_luna = 10 ** (-0.4 * (3.84 + 0.026 * abs(fase_deg) + 4e-9 * fase_deg ** 4))
        f_rho = 10 ** 5.36 * (1.06 + math.cos(rho * D2R) ** 2) + 10 ** (6.15 - rho / 40)
        luz_v = f_rho * i_luna * 10 ** (-0.4 * kv * x_ks(zm)) * (1 - 10 ** (-0.4 * kv * x_ks(z)))
    res = {}
    for filt in "BV":
        k = te.EXTINCION[filt]
        oscuro = CIELO_OSCURO[filt] - 2.5 * math.log10(x_ks(z)) + k * (x_ks(z) - 1)
        flujo = 10 ** (-0.4 * oscuro)
        if luz_v > 0:
            flujo += 10 ** (-0.4 * (nl_a_mag(luz_v) + (COLOR_LUZ_LUNA if filt == "B" else 0.0)))
        res[filt] = -2.5 * math.log10(flujo)
    return res


def estado(nombre, fecha, hora):
    jd = jd_local(fecha, hora)
    lam, beta, paralaje, T = luna(jd)
    ra_l, de_l = ecl_a_ecuatorial(lam, beta, T)
    alt_l = altura(ra_l, de_l, jd)
    alt_l -= paralaje * math.cos(alt_l * D2R)                    # altura topocéntrica
    elong = math.acos(math.cos(beta * D2R) * math.cos((lam - sol(jd)) * D2R)) / D2R
    fase = 180 - elong                                          # ángulo de fase (0 = llena)
    ilum = (1 + math.cos(fase * D2R)) / 2
    ra, dec = OBJETOS[nombre]
    alt = altura(ra, dec, jd)
    rho = separacion(ra, dec, ra_l, de_l)
    return {"alt": alt, "alt_luna": alt_l, "ilum": ilum, "rho": rho, "fase": fase,
            "ra_luna": ra_l, "dec_luna": de_l,
            "cielo": brillo_cielo(alt, alt_l, rho, fase) if alt > 5 else None}


def snr_objeto(nombre, alt, cielo, inst=te.INSTRUMENTO):
    """S/N al brillo superficial medio, con el plan y el método de tiempos_exposicion.py."""
    v, b, (a, c) = next((v, b, tam) for n, v, b, tam, _ in te.OBJETOS if n == nombre)
    tb, nb, tv, nv = te.PLAN.get(nombre, te.PLAN_DEFECTO)
    x = te.airmass(alt)
    area_obj = math.pi / 4 * a * c * 3600
    area_res = math.pi / 4 * te.DIAM_SN ** 2
    npix = area_res / inst["escala"] ** 2
    out = {}
    for filt, m, t, n in (("B", b, tb, nb), ("V", v, tv, nv)):
        zp = te.zp_efectivo(filt, x, inst)
        mu = m + 2.5 * math.log10(area_obj)
        cielo_px = te.tasa(cielo[filt], zp) * inst["escala"] ** 2
        out[filt] = (te.snr(te.tasa(mu, zp) * area_res, cielo_px, npix, t, n, inst), cielo_px * t)
    return out


def main():
    for noche, fecha in NOCHES.items():
        e0 = estado("NGC 7793", fecha, 22)
        print(f"\n## Noche del {noche}: Luna al {100 * e0['ilum']:.0f} % "
              f"(a las 22:00 en AR {e0['ra_luna'] / 15:.2f} h, Dec {e0['dec_luna']:+.1f}°)\n")
        print("| Hora | Alt. Luna | Objeto | Alt. | Dist. Luna | μB | μV | S/N B | S/N V | Cielo B (e⁻/px por exp.) |")
        print("|---|---|---|---|---|---|---|---|---|---|")
        for hora in HORAS:
            for nombre in OBJETOS:
                e = estado(nombre, fecha, hora)
                if e["alt"] < 30:
                    continue
                s = snr_objeto(nombre, e["alt"], e["cielo"])
                print(f"| {hora % 24:02d}:00 | {e['alt_luna']:.0f}° | {nombre} | {e['alt']:.0f}° | "
                      f"{e['rho']:.0f}° | {e['cielo']['B']:.1f} | {e['cielo']['V']:.1f} | "
                      f"{s['B'][0]:.0f} | {s['V'][0]:.0f} | {s['B'][1]:.0f} |")
    print("\nReferencia sin Luna (8/10, cielo de tiempos_exposicion.py): NGC 7793 S/N B 20 · V 22; "
          "NGC 1291 S/N B 13 · V 19.")


if __name__ == "__main__":
    main()
