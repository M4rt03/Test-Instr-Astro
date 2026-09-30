# Candidatos para la Tarea 2 (2026): objetos y estrellas estándar

Listas para elegir el objeto de cada grupo y las estrellas estándar comunes
para las noches de observación del **jueves 8 y el jueves 15 de octubre de 2026**
con el telescopio MAS500 (campo de 24,6′ × 24,6′, filtros B y V).

## Criterios del enunciado

- El objeto debe ser una **nebulosa o galaxia**, con magnitud aparente entre **~6 y 10**
  en la banda correspondiente.
- Debe ser **visible ambas noches antes de la 01:00** desde Cerro Tololo.
- Las estándares deben tener magnitud **entre 6 y 8 en cada banda** (B y V).

## Cómo se calculó

- **Sitio:** Cerro Tololo (−30,169°, −70,807°, 2200 m). La hora local es UTC−3.
- **Ventana útil:** desde el fin del crepúsculo náutico (Sol a −12°, ~20:50) hasta la 01:00.
- **"h > 30°"** son las horas de esa ventana con altura mayor a 30° (masa de aire < 2).
- **Fuentes:**
  - Objetos: [OpenNGC](https://github.com/mattiaverga/OpenNGC) y RC3 (VizieR VII/155).
  - Estrellas: Hipparcos (VizieR I/239), con magnitudes y tipos espectrales de SIMBAD.
- **Verificación de las estrellas:** son estrellas simples en SIMBAD, sin bandera de
  variabilidad (`VarFlag`) ni solución de componentes (`MultFlag`) en Hipparcos.

## Condiciones de las noches

| Noche | Luna | Consecuencia |
|---|---|---|
| Jue 8/10/2026 | Nueva (3 % iluminada), bajo el horizonte | Cielo oscuro toda la ventana |
| Jue 15/10/2026 | Creciente (26 %), en Sagitario, se pone ~00:50 | Evitar objetos en Sagitario: queda a ~15° de M17 y a ~7° de M8 y M20 |

## 1. Objetos candidatos (10)

| # | Objeto | Tipo | AR (J2000) | Dec (J2000) | B | V | Tamaño | Alt. máx | h > 30° 8/10 · 15/10 | Dist. Luna 15/10 | Notas |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **NGC 253** (galaxia del Escultor) | galaxia espiral | 00:47:33 | −25:17:17 | 8,0 | 7,2 | 27′ × 5′ | 83° | 4,3 · 4,3 | 93° | Brillante y con mucha estructura (polvo, brotes de formación estelar). Estándar cercana: HIP 4442 |
| 2 | **NGC 7293** (nebulosa de la Hélice) | nebulosa planetaria | 22:29:38 | −20:50:14 | 7,5 | 7,3 | 16′ | 80° | 4,3 · 4,3 | 67° | Cabe completa en el campo. Brillo superficial bajo: requiere exposiciones largas |
| 3 | **NGC 55** | galaxia de canto | 00:14:53 | −39:11:47 | 8,6 | 7,9 | 30′ × 3′ | 81° | 4,3 · 4,3 | 80° | Más larga que el lado del campo: cabe orientada en la diagonal |
| 4 | **NGC 300** | galaxia espiral | 00:54:53 | −37:41:03 | 8,7 | 8,1 | 19′ × 13′ | 81° | 4,3 · 4,3 | 88° | Brillo superficial bajo |
| 5 | **NGC 7793** | galaxia espiral | 23:57:49 | −32:35:27 | 9,6 | 9,1 | 10′ × 6′ | 88° | 4,3 · 4,3 | 80° | Pasa casi por el cenit; tamaño cómodo para el campo |
| 6 | **NGC 7009** (nebulosa de Saturno) | nebulosa planetaria | 21:04:10 | −11:21:47 | 8,3 | 8,0 | 0,7′ | 71° | 4,3 · 4,3 | 52° | Pequeña y de alto brillo superficial; el campo tiene muchas estrellas para medir el seeing |
| 7 | **NGC 6744** | galaxia espiral | 19:09:46 | −63:51:27 | 9,1 | 8,3 | 16′ × 10′ | 55° | 4,3 · 4,0 | 40° | Circumpolar; mejor al comienzo de la noche |
| 8 | **NGC 6818** (Little Gem) | nebulosa planetaria | 19:43:57 | −14:09:11 | 9,9 | 9,3 | 0,8′ | 72° | 3,8 · 3,3 | 33° | Observar temprano (va bajando durante la noche) |
| 9 | **NGC 246** (nebulosa de la Calavera) | nebulosa planetaria | 00:47:03 | −11:52:19 | 8,0 | 10,9 | 4′ | 71° | 3,8 · 4,3 | 100° | V en el límite del rango; comparte estándar con NGC 253 |
| 10 | **NGC 6618** (M17, nebulosa Omega) | región HII | 18:20:47 | −16:10:17 | 6,0 | 7,0 | 13′ | 60° | 2,5 · 2,0 | 15° | Solo al comienzo de la noche, y la Luna molesta el 15/10. Es el objeto del informe de 2025 |

**Reemplazos posibles:**
- NGC 247 (galaxia; B 9,7, V 9,2; 20′ × 6′).
- NGC 6822 (galaxia de Barnard; B 9,4, V 10,0; 17′).
- NGC 1097 (galaxia; B 10,1, V 9,8; sale tarde, a 60° de altura a la 01:00).

**No recomendados:**
- M8 (Laguna) y M20 (Trífida): quedan a ~7° de la Luna el 15/10, y M8 es más grande que el campo.

### Origen de las magnitudes

- **NGC 253 y NGC 7793:** RC3, con V = B_T − (B−V)_T. OpenNGC da V = 11,1 para NGC 253, que es un error.
- **NGC 55, NGC 300, NGC 6744 y la V de NGC 6818:** SIMBAD.
- **Nebulosas planetarias y M17:** magnitudes **integradas** de OpenNGC. SIMBAD les asigna la
  magnitud de la estrella central (B ≈ 11,5–13,5), que **no** es la de la nebulosa. Hay que
  aclarar en el informe cuál se usa.

## 2. Estrellas estándar candidatas (8)

Son enanas de tipo A0–A9 cerca del cenit (δ ≈ −22° a −31°), repartidas en ascensión recta
para que siempre haya una alta a cualquier hora de la ventana.

| # | Estrella | HD | AR (J2000) | Dec (J2000) | V | B | B−V | Tipo | Alt. máx | h > 30° 8/10 · 15/10 | Sirve para |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **HIP 97210** | 186417 | 19:45:24 | −30:54:11 | 6,81 | 6,93 | 0,12 | A1V | 83° | 4,3 · 3,8 | Comienzo de noche: NGC 6744, NGC 6818, M17 |
| 2 | **HIP 98926** | 190285 | 20:05:12 | −26:48:48 | 7,22 | 7,31 | 0,09 | A1V | 86° | 4,3 · 4,0 | Comienzo de noche |
| 3 | **HIP 105318** | 202941 | 21:19:51 | −27:12:33 | 7,05 | 7,05 | 0,00 | A0V | 87° | 4,3 · 4,3 | NGC 7009 |
| 4 | **HIP 107766** | 207480 | 21:49:54 | −27:24:14 | 7,11 | 7,17 | 0,06 | A1V | 87° | 4,3 · 4,3 | NGC 7009, NGC 7293 |
| 5 | **HIP 110861** | 212852 | 22:27:39 | −26:25:08 | 7,20 | 7,41 | 0,21 | A9V | 86° | 4,3 · 4,3 | NGC 7293 (a 6°) |
| 6 | **HIP 116750** | 222332 | 23:39:42 | −22:32:00 | 7,29 | 7,41 | 0,12 | A0V | 82° | 4,3 · 4,3 | NGC 7793, NGC 55 |
| 7 | **HIP 4442** | 5524 | 00:56:49 | −25:21:48 | 7,24 | 7,33 | 0,09 | A2/3V | 82° | 4,2 · 4,3 | NGC 253 (a 2°), NGC 246, NGC 300 |
| 8 | **HIP 6897** | 9063 | 01:28:48 | −24:47:53 | 7,04 | 7,27 | 0,23 | A8V | 75° | 3,7 · 4,0 | Final de la noche |

Las magnitudes B y V son las de SIMBAD. Las 8 cumplen 6 ≤ B, V ≤ 8.

**Descartadas en la verificación:**
- **HIP 112362 y HIP 6393:** Hipparcos las marca como posibles variables (`VarFlag = 1`).
- **HIP 4496 y HIP 6507:** Hipparcos las resuelve como dobles (`MultFlag = C`).
- **HIP 116375 y HIP 117678** (las del informe de 2025): son gigantes K con B ≈ 8,1, fuera del
  rango 6–8 en B.

## Recomendaciones para el plan de observación

- **Estándar temprana + estándar tardía:** elegir una del grupo 1–2 y otra del 7–8. Así se
  observan a masas de aire distintas y se puede **medir** el coeficiente de extinción k.
  Si no, el pipeline usa el valor típico de Cerro Tololo (`extincion` en `config.yaml`).
- **Tiempos de exposición de las estándares:** con V ≈ 7, conviene partir las pruebas en 1–3 s
  y revisar que el pico de la estrella quede bajo ~60 000 ADU. El pipeline excluye del
  zeropoint las mediciones saturadas.
- **Verificación en SIMBAD:** antes de fijar el plan, conviene revisar la ficha de cada
  estrella elegida (https://simbad.cds.unistra.fr) por si hay información nueva.
