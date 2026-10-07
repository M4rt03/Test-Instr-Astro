"""Mide la escala de placa del MAS500 comparando estrellas detectadas con Gaia DR3.

Usa las 8 fuentes detectadas en una exposición de 9 s en B de la estrella HIP 117445
(2023-10-06, tabla impresa en Calibracao.ipynb). Descarga de VizieR las estrellas de
Gaia DR3 del campo y busca la escala, la rotación y la orientación que hacen coincidir
las posiciones. Luego ajusta una transformación lineal por mínimos cuadrados.

Solo usa la biblioteca estándar:

    python medir_escala.py
"""

import math
import urllib.request

D2R = math.pi / 180

# (x, y) de SourceCatalog en Calibracao.ipynb; la fuente 5 es HIP 117445
DETECCIONES = [
    (3166.708, 14.210), (3777.846, 1011.308), (1939.996, 1095.446), (2731.252, 1607.438),
    (2042.088, 2036.037), (3709.030, 2703.857), (467.218, 3434.992), (3191.740, 4092.187),
]
ANCLA = 4
HIP117445 = (357.20569, 2.21441)            # SIMBAD (HD 223346), J2000
EPOCA = 2023.76                              # DATE-OBS = 2023-10-06T03:32:46 UTC
PIXEL_UM, FOCAL_MM = 9.0, 6500.0             # XPIXSZ y FOCALLEN del header


def gaia(ra, dec, radio_arcmin=19, gmax=15):
    url = ("https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=I/355/gaiadr3"
           f"&-c={ra}{dec:+}&-c.rm={radio_arcmin}&-out=RA_ICRS,DE_ICRS,pmRA,pmDE,Gmag"
           f"&Gmag=%3C{gmax}&-out.max=2000")
    lineas = urllib.request.urlopen(url, timeout=120).read().decode().splitlines()
    datos = [l for l in lineas if l and not l.startswith("#")][3:]   # salta encabezados
    estrellas = []
    for l in datos:
        p = l.split("\t")
        ra_, de_ = float(p[0]), float(p[1])
        pmra = float(p[2]) if p[2].strip() else 0.0
        pmde = float(p[3]) if p[3].strip() else 0.0
        dt = EPOCA - 2016.0                          # época de Gaia DR3
        estrellas.append((ra_ + pmra / 3.6e6 * dt / math.cos(de_ * D2R), de_ + pmde / 3.6e6 * dt))
    return estrellas


def plano_tangente(ra, de, ra0, de0):
    """Proyección gnomónica en arcsec alrededor de (ra0, de0)."""
    a, d, a0, d0 = ra * D2R, de * D2R, ra0 * D2R, de0 * D2R
    c = math.sin(d0) * math.sin(d) + math.cos(d0) * math.cos(d) * math.cos(a - a0)
    xi = math.cos(d) * math.sin(a - a0) / c
    eta = (math.cos(d0) * math.sin(d) - math.sin(d0) * math.cos(d) * math.cos(a - a0)) / c
    return xi / D2R * 3600, eta / D2R * 3600


def coincidencias(off, cat, s, th, par, tol=6.0):
    ct, st = math.cos(th), math.sin(th)
    n, err, pares = 0, 0.0, []
    for dx, dy in off:
        u = par * dx
        xi, eta = s * (ct * u - st * dy), s * (st * u + ct * dy)
        j = min(range(len(cat)), key=lambda k: math.hypot(xi - cat[k][0], eta - cat[k][1]))
        d = math.hypot(xi - cat[j][0], eta - cat[j][1])
        if d < tol:
            n, err = n + 1, err + d
            pares.append(j)
        else:
            pares.append(None)
    return n, err, pares


def minimos_cuadrados(filas, b):
    """Resuelve b = A·p (A de n×3) por ecuaciones normales."""
    ata = [[sum(f[i] * f[j] for f in filas) for j in range(3)] for i in range(3)]
    atb = [sum(f[i] * v for f, v in zip(filas, b)) for i in range(3)]
    m = [ata[i] + [atb[i]] for i in range(3)]
    for i in range(3):
        piv = max(range(i, 3), key=lambda r: abs(m[r][i]))
        m[i], m[piv] = m[piv], m[i]
        for r in range(3):
            if r != i:
                fac = m[r][i] / m[i][i]
                m[r] = [m[r][k] - fac * m[i][k] for k in range(4)]
    return [m[i][3] / m[i][i] for i in range(3)]


def main():
    cat = [plano_tangente(r, d, *HIP117445) for r, d in gaia(*HIP117445)]
    ax, ay = DETECCIONES[ANCLA]
    otras = [p for i, p in enumerate(DETECCIONES) if i != ANCLA]
    off = [(x - ax, y - ay) for x, y in otras]

    # 1) Búsqueda en grilla: escala 0.24-0.40 "/px, rotación 0-360°, con y sin reflexión
    mejor = max(((*coincidencias(off, cat, s / 10000, t / 10 * D2R, par)[:2], s / 10000, t / 10, par)
                 for par in (1, -1) for s in range(2400, 4001, 4) for t in range(0, 3600, 5)),
                key=lambda r: (r[0], -r[1]))
    n, _, s, th, par = mejor
    print(f"Búsqueda: {n} de {len(off)} estrellas coinciden con s = {s:.4f}\"/px, rotación {th:.1f}°")
    for s_prueba in (206265 * PIXEL_UM * 1e-3 / FOCAL_MM, 0.36):
        n_p = max(coincidencias(off, cat, s_prueba + ds / 10000, t / 10 * D2R, p)[0]
                  for p in (1, -1) for ds in range(-15, 16, 3) for t in range(0, 3600, 5))
        print(f"  con {s_prueba:.4f}\"/px coinciden como máximo {n_p} de {len(off)}")

    # 2) Ajuste lineal con todas las coincidencias (incluida HIP 117445 en el origen)
    _, _, pares = coincidencias(off, cat, s, th * D2R, par)
    filas, xis, etas = [[ax, ay, 1.0]], [0.0], [0.0]
    for (x, y), j in zip(otras, pares):
        if j is not None:
            filas.append([x, y, 1.0])
            xis.append(cat[j][0])
            etas.append(cat[j][1])
    a, b, c = minimos_cuadrados(filas, xis)
    d, e, f = minimos_cuadrados(filas, etas)
    escala = math.sqrt(abs(a * e - b * d))
    res = [math.hypot(a * x + b * y + c - xi, d * x + e * y + f - eta)
           for (x, y, _), xi, eta in zip(filas, xis, etas)]
    print(f"\nAjuste con {len(filas)} estrellas:")
    print(f"  escala      = {escala:.4f}\"/px  (x: {math.hypot(a, d):.4f}, y: {math.hypot(b, e):.4f})")
    print(f"  rotación    = {math.degrees(math.atan2(d, a)):.2f}°")
    print(f"  residuo rms = {math.sqrt(sum(r * r for r in res) / len(res)):.2f}\"")
    print(f"  campo       = {4096 * escala / 60:.2f}' × {4096 * escala / 60:.2f}'")
    print(f"  focal efectiva = {206265 * PIXEL_UM * 1e-3 / escala:.0f} mm "
          f"(header: {FOCAL_MM:.0f} mm, que daría {206265 * PIXEL_UM * 1e-3 / FOCAL_MM:.4f}\"/px)")


if __name__ == "__main__":
    main()
