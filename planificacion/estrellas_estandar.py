"""Magnitudes Johnson B y V de las estrellas estándar, consultadas con astroquery en VizieR.

Para cada estrella de candidatos_2026.md (8 seleccionadas, 3 de reserva y
las 2 reemplazadas por tener B−V de Tycho) consulta:

- II/168/ubvmeans: Mermilliod (1991), medias homogéneas de fotometría fotoeléctrica en el
  sistema UBV de Johnson. Es la referencia para estrellas brillantes.
- I/239/hip_main: Hipparcos. Trae el origen de cada valor: r_Vmag y r_B-V = G (fotometría
  terrestre, Johnson), H (V calculada desde Hp) o T (calculada desde Tycho).

Elige en este orden: Mermilliod; si no está, Hipparcos. La calidad queda "baja" cuando el
B−V viene de Tycho: es una transformación, no una medición en Johnson. APASS (II/336) no se
usa porque satura en estrellas de V < ~10 (ver calibraciones.md).

Imprime una tabla con ambas fuentes, los tiempos de exposición (el mismo modelo de
tiempos_exposicion.py) y el bloque de estándares para pipeline/config.yaml.

Requiere astroquery (pip install astroquery):

    python estrellas_estandar.py
    python estrellas_estandar.py --todas --csv estandares_bv.csv
"""

import argparse
import csv
import math

from astroquery.vizier import Vizier

import tiempos_exposicion as te

SEL, RES, EXT = "seleccionada", "reserva", "solo extinción"
ESTRELLAS = [  # (HD, uso)
    (202941, SEL), (207480, SEL), (210300, SEL), (212643, SEL),
    (220881, SEL), (562, SEL), (8130, SEL), (12206, SEL),
    (223884, RES), (225200, RES), (7323, RES),
    (195500, EXT), (215863, EXT),         # B−V de Tycho: solo para medir la extinción
]
DIFERENCIA_MAX = 0.03   # mag: avisar si Mermilliod e Hipparcos difieren más que esto


def _num(v):
    """Valor numérico de una celda de astropy, o None si está vacía (enmascarada)."""
    if v is None or bool(getattr(v, "mask", False)):
        return None
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(x) else x


def _txt(v):
    if v is None or bool(getattr(v, "mask", False)):
        return ""
    return str(v).strip()


def consultar(hd):
    """Filas de Hipparcos y de Mermilliod para HD `hd` (None si no aparece)."""
    hip = Vizier(columns=["HD", "HIP", "Vmag", "r_Vmag", "B-V", "e_B-V", "r_B-V",
                          "HvarType", "SpType", "RAICRS", "DEICRS"], row_limit=-1)
    t = hip.query_constraints(catalog="I/239/hip_main", HD=str(hd))
    fila_hip = t[0][0] if len(t) and len(t[0]) else None

    # En Mermilliod, el identificador de una estrella HD es "+100" seguido del número
    # con 6 dígitos (HD 562 → +100000562). Se prefiere la fila sin componente (m_LID vacío).
    ubv = Vizier(columns=["LID", "m_LID", "Vmag", "e_Vmag", "o_Vmag", "B-V", "e_B-V", "o_B-V"],
                 row_limit=-1)
    t = ubv.query_constraints(catalog="II/168/ubvmeans", LID=f"+100{hd:06d}")
    fila_ubv = None
    if len(t) and len(t[0]):
        filas = sorted(t[0], key=lambda f: _txt(f["m_LID"]) != "")
        fila_ubv = filas[0]
    return fila_hip, fila_ubv


def elegir(hd, fila_hip, fila_ubv):
    r = {"nombre": f"HD {hd}", "HIP": None, "ra": None, "dec": None, "tipo": "",
         "V_merm": None, "BV_merm": None, "n_merm": None,
         "V_hip": None, "BV_hip": None, "origen_hip": ""}
    if fila_hip is not None:
        r.update(HIP=int(fila_hip["HIP"]), ra=_num(fila_hip["RAICRS"]), dec=_num(fila_hip["DEICRS"]),
                 tipo=_txt(fila_hip["SpType"]), V_hip=_num(fila_hip["Vmag"]),
                 BV_hip=_num(fila_hip["B-V"]),
                 origen_hip=f"{_txt(fila_hip['r_Vmag'])}/{_txt(fila_hip['r_B-V'])}")
    if fila_ubv is not None:
        r.update(V_merm=_num(fila_ubv["Vmag"]), BV_merm=_num(fila_ubv["B-V"]),
                 n_merm=int(_num(fila_ubv["o_Vmag"]) or 0))

    if r["V_merm"] is not None and r["BV_merm"] is not None:
        v, bv = r["V_merm"], r["BV_merm"]
        r["fuente"] = f"Mermilliod ({r['n_merm']} obs.)"
        r["calidad"] = "alta" if r["n_merm"] >= 3 else "media"
    elif r["V_hip"] is not None and r["BV_hip"] is not None:
        v, bv = r["V_hip"], r["BV_hip"]
        r["fuente"] = f"Hipparcos ({r['origen_hip']})"
        if r["origen_hip"].endswith("T"):
            r["calidad"] = "baja"                 # B−V transformado desde Tycho
        elif r["origen_hip"] == "G/G":
            r["calidad"] = "alta"
        else:
            r["calidad"] = "media"                # V calculada desde Hp
    else:
        raise RuntimeError(f"HD {hd}: no hay V y B−V en Mermilliod ni en Hipparcos")
    r["V"], r["B"] = round(v, 3), round(v + bv, 3)

    r["aviso"] = ""
    if None not in (r["V_merm"], r["V_hip"], r["BV_merm"], r["BV_hip"]):
        dv, dbv = r["V_merm"] - r["V_hip"], r["BV_merm"] - r["BV_hip"]
        if max(abs(dv), abs(dbv)) > DIFERENCIA_MAX:
            r["aviso"] = f"Mermilliod − Hipparcos: ΔV = {dv:+.3f}, Δ(B−V) = {dbv:+.3f}"
    return r


def sexagesimal(grados, horas=False):
    x = grados / 15 if horas else grados
    signo = "-" if x < 0 else ("" if horas else "+")
    decimas = round(abs(x) * 36000)          # en décimas de segundo, para no imprimir 60.0
    if horas:
        decimas %= 24 * 36000
    a, resto = divmod(decimas, 36000)
    b, c = divmod(resto, 600)
    return f"{signo}{a:02d}:{b:02d}:{c / 10:04.1f}"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--todas", action="store_true", help="incluir las de reserva en el bloque YAML")
    p.add_argument("--csv", help="guardar la tabla en este archivo CSV")
    p.add_argument("--noche", default="noche1", help="carpeta de la noche para las rutas del YAML")
    args = p.parse_args()

    inst = te.INSTRUMENTO
    estrellas = []
    for hd, uso in ESTRELLAS:
        r = elegir(hd, *consultar(hd))
        r["uso"] = uso
        for filt in "BV":
            r[f"t_{filt}"] = te.t_estandar(r[filt], filt, inst)[0]
        estrellas.append(r)

    print("| Estrella | Uso | Tipo | V Merm. | B−V Merm. | V Hip. | B−V Hip. (origen) | "
          "Fuente elegida | Calidad | B | V | t_exp B | t_exp V |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    fmt = lambda x: "—" if x is None else f"{x:.3f}"
    for r in estrellas:
        print(f"| {r['nombre']} | {r['uso']} | {r['tipo']} | "
              f"{fmt(r['V_merm'])} | {fmt(r['BV_merm'])} | {fmt(r['V_hip'])} | "
              f"{fmt(r['BV_hip'])} ({r['origen_hip']}) | {r['fuente']} | {r['calidad']} | "
              f"{r['B']:.3f} | {r['V']:.3f} | {r['t_B']:.1f} s | {r['t_V']:.1f} s |")
    avisos = [r for r in estrellas if r["aviso"]]
    for r in avisos:
        print(f"\nAVISO {r['nombre']}: {r['aviso']}")
    print("\nOrigen en Hipparcos (V/B−V): G = fotometría terrestre, H = desde Hp, T = desde Tycho.")

    print(f"\n# Bloque para pipeline/config.yaml (dentro de la noche, en 'estandares:')")
    for r in estrellas:
        if not (r["uso"] == SEL or args.todas):
            continue
        carpeta = r["nombre"].replace(" ", "_")
        print(f'      - nombre: "{r["nombre"]}"\n'
              f'        ra: "{sexagesimal(r["ra"], horas=True)}"\n'
              f'        dec: "{sexagesimal(r["dec"])}"\n'
              f'        magnitudes: {{B: {r["B"]:.3f}, V: {r["V"]:.3f}}}   # {r["fuente"]}, calidad {r["calidad"]}\n'
              f'        archivos:\n'
              f'          B: "datos/{args.noche}/{carpeta}/*_B_*.fit"\n'
              f'          V: "datos/{args.noche}/{carpeta}/*_V_*.fit"')

    if args.csv:
        campos = ["nombre", "HIP", "uso", "tipo", "ra", "dec", "V_merm", "BV_merm", "n_merm",
                  "V_hip", "BV_hip", "origen_hip", "fuente", "calidad", "B", "V", "t_B", "t_V", "aviso"]
        with open(args.csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=campos, extrasaction="ignore")
            w.writeheader()
            w.writerows(estrellas)
        print(f"\nTabla guardada en {args.csv}")


if __name__ == "__main__":
    main()
