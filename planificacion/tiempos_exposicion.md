# Tiempos de exposición para el MAS500 (Tarea 2, 2026)

Estimación previa de los tiempos de exposición para los objetos y las estrellas
estándar de [`candidatos_2026.md`](candidatos_2026.md). Las tablas se generan con
[`tiempos_exposicion.py`](tiempos_exposicion.py), que solo usa la biblioteca estándar de Python:

```bash
python tiempos_exposicion.py
python tiempos_exposicion.py --zp-v 21.8 --ruido 3.2 --seeing 3.5   # con valores medidos
```

## Qué se pide

"Calcular el tiempo de exposición" significa decidir, para cada objeto y cada filtro:

1. **Cuánto dura cada exposición (t_exp).** Tiene un máximo: ni el objeto ni las estrellas
   que se usarán deben saturar. Tiene también un mínimo: si la exposición es muy corta, domina
   el ruido de lectura del detector.
2. **Cuántas exposiciones tomar (N).** Con eso se alcanza la señal a ruido (S/N) buscada en el
   apilado, dentro del bloque de 30 minutos del plan.

Es lo que pide el punto 7 del enunciado: llegar con una estimación, confirmarla con
exposiciones de prueba al comienzo de cada bloque y **justificar en el informe** la elección final.

## Cómo se calcula

### Ecuación de señal a ruido

La S/N de una fuente medida en `n_pix` píxeles, apilando N exposiciones de t segundos, es:

```
              N · S · t
S/N = ─────────────────────────────────────────
      √( N · [ S·t + n_pix·(C·t + D·t + RN²) ] )
```

| Símbolo | Significado | Unidad |
|---|---|---|
| S | señal de la fuente | e⁻/s |
| C | cielo por píxel | e⁻/s/px |
| D | corriente oscura | e⁻/s/px |
| RN | ruido de lectura | e⁻ |

Las tasas se obtienen de las magnitudes con el zeropoint (ZP, la magnitud que da 1 e⁻/s) y la
extinción: `S = 10^(0.4·(ZP − k·X − m))`. Para objetos extendidos se usa el brillo superficial
medio, `μ = m + 2.5 log10(área en arcsec²)`, y la S/N se calcula en un elemento de seeing
(un círculo de diámetro igual al FWHM).

### Saturación

- **Estrella:** el píxel central recibe una fracción `1/(1.133·FWHM_px²)` del flujo total.
  Satura cuando esa fracción, sumada al cielo, supera el límite de linealidad.
- **Objeto extendido:** un píxel satura cuando `10^(0.4·(ZP − μ))·escala²·t` llega a ese límite.

### Parámetros usados

| Parámetro | Valor | Origen |
|---|---|---|
| Escala de placa | 0,317″/px | Medida con estrellas de Gaia ([`medir_escala.py`](medir_escala.py): 0,3168″/px); coincide con astrometry.net en 2025 (0,3163″/px). El 0,36″/px que usó el informe de 2025 era incorrecto |
| Ganancia | 0,85 e⁻/ADU | Header `GAIN` de los datos de 2025 |
| Ruido de lectura | 3,9 e⁻ | Ficha de Moravian: GSENSE4040 FSI en modo 16-bit HDR (el de ganancia 0,85 e⁻/ADU). Conviene verificarlo en los bias |
| Corriente oscura | 0,05 e⁻/s/px | Supuesto: Moravian no lo publica. **Medir** en los darks |
| Límite de saturación | 40 000 ADU sobre el bias | Conservador. El ADC satura en 65 535 ADU, antes que el pozo (56 600 e⁻ = 66 600 ADU); según Moravian, el modo HDR es lineal en todo el rango |
| ZP sobre la atmósfera | B 21,98 · V 22,05 (e⁻/s) | Estimado; ver detalle debajo |
| Extinción | k_B = 0,25 · k_V = 0,15 | Valores típicos de Cerro Tololo |
| Cielo sin Luna (8/10) | B 22,7 · V 21,8 mag/arcsec² | Valores típicos de Cerro Tololo |
| Cielo con Luna al 25 % (15/10) | B 22,1 · V 21,6 mag/arcsec² | Valores típicos para esa fase |
| Seeing | 4,0″ para la S/N · 2,5″ para la saturación | El informe 2025 midió 3,9–4,4″; la saturación se evalúa con buen seeing |

**Cómo se estimó el ZP:**
- área efectiva de 1570 cm² (espejo de 0,5 m con una obstrucción de ~45 %);
- flujo de Vega en B y V (Bessell et al. 1998);
- eficiencia total de 0,30 en B y 0,48 en V (espejos, filtro y eficiencia cuántica).

**El ZP es la mayor incertidumbre: puede estar errado en ±0,5 mag.** Con un ZP 0,5 mag peor, la
S/N de los objetos débiles baja en un factor ~1,5 (NGC 300 en B: 13,6 → 9,0).

### Cómo medir los valores reales con los datos de 2025

Con el notebook `pipeline/notebooks/tarea2_paso_a_paso.ipynb`, ejecutando solo
`paso_calibracion` y `paso_estandares` (no necesitan las imágenes de ciencia), se obtiene:

- **`calibracion.json`:** `ruido_lectura_e`, `dark_current_mediana_adu_s` (multiplicar por la
  ganancia para pasar a e⁻/s) y `bias_mediana_adu`.
- **`zeropoints.ecsv`:** el ZP en ADU/s, ya corregido por extinción. Para el script, el ZP en
  e⁻/s es `ZP_ADU + 2.5·log10(0,85)` ≈ `ZP_ADU − 0,18`.
- **`estandares.ecsv`:** el FWHM de cada exposición, para `--seeing`.

Luego: `python tiempos_exposicion.py --zp-b … --zp-v … --ruido … --oscuridad … --seeing …`

## Resultado principal: dominan el ruido de lectura y la corriente oscura

Con píxeles de 0,317″, cada píxel recibe muy poca luz del cielo: entre 0,04 y 0,13 e⁻/s.
Para que el ruido del cielo domine sobre el ruido de lectura (cielo > 10·RN² por píxel) harían
falta exposiciones de 20 a 65 minutos. Por lo tanto:

- **Conviene tomar pocas exposiciones largas en lugar de muchas cortas.** Con el mismo
  tiempo total (~14 min en B, NGC 300):

  | Exposición | S/N |
  |---|---|
  | 30 s × 28 | 8,4 |
  | 60 s × 14 | 11,0 |
  | 120 s × 7 | 13,6 |
  | 300 s × 3 | 17,1 |

- **El límite práctico de t_exp lo ponen otras cosas:** el seguimiento o guiado de la montura,
  los satélites (con pocas exposiciones, perder una cuesta más) y la saturación del núcleo del
  objeto. Por eso se propone 90–120 s y no 300 s. Si el guiado funciona bien, se puede subir.
- **La Luna del 15/10 casi no cambia la S/N** (NGC 300 en B: 13,6 → 12,8), porque el cielo no
  es la fuente dominante de ruido. Sí cambia el nivel de fondo, que debe sustraerse bien.
- **Los 30 s de 2025 (M17) quedaban en el régimen dominado por el ruido de lectura.** Para M17 no
  importaba, porque es muy brillante. Para galaxias de ~23 mag/arcsec² sí importa.

## 1. Objetos: plan por bloque de 30 min

Los tiempos suman ~24 min de exposición más lecturas, lo que deja ~6 min para apuntar, enfocar
y tomar las pruebas. B recibe más tiempo porque el cielo es más oscuro y la eficiencia, menor.

| Objeto | X | B: t_exp × N | V: t_exp × N | Tiempo total (min) | μB medio | μV medio | S/N B 8/10 · 15/10 | S/N V 8/10 · 15/10 | S/N B borde (μ+2) 8/10 | μ límite B 8/10 (S/N=3) | μ satura B · V | Estrella satura V |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| NGC 7293 | 1.04 | 120 s × 7 | 90 s × 6 | 24 | 22.2 | 22.0 | 40 · 38 | 37 · 36 | 7 | 25.1 | 13.1 · 13.0 | V < 10.8 |
| NGC 7009 | 1.06 | 30 s × 12 | 20 s × 12 | 12 | 16.2 | 15.9 | 855 · 855 | 870 · 870 | 313 | 24.1 | 11.6 · 11.3 | V < 9.2 |
| NGC 300 | 1.06 | 120 s × 7 | 90 s × 6 | 24 | 23.4 | 23.3 | 14 · 13 | 12 · 11 | 2 | 25.1 | 13.1 · 13.0 | V < 10.8 |
| NGC 7793 | 1.04 | 120 s × 7 | 90 s × 6 | 24 | 22.8 | 22.4 | 24 · 22 | 26 · 25 | 4 | 25.1 | 13.1 · 13.0 | V < 10.8 |
| NGC 247 | 1.06 | 120 s × 7 | 90 s × 6 | 24 | 23.4 | 22.9 | 14 · 13 | 17 · 16 | 2 | 25.1 | 13.1 · 13.0 | V < 10.8 |
| NGC 1097 | 1.15 | 90 s × 9 | 60 s × 9 | 24 | 23.3 | 23.0 | 14 · 13 | 14 · 13 | 2 | 25.0 | 12.8 · 12.5 | V < 10.4 |
| NGC 1291 | 1.22 | 90 s × 9 | 60 s × 9 | 24 | 23.1 | 22.4 | 16 · 15 | 22 · 22 | 3 | 25.0 | 12.7 · 12.5 | V < 10.4 |
| NGC 1316 | 1.22 | 90 s × 9 | 60 s × 9 | 24 | 23.1 | 22.2 | 17 · 16 | 28 · 28 | 3 | 25.0 | 12.7 · 12.5 | V < 10.4 |
| NGC 1313 | 1.41 | 120 s × 7 | 90 s × 6 | 24 | 23.3 | 23.1 | 14 · 13 | 13 · 13 | 2 | 25.0 | 13.0 · 12.9 | V < 10.8 |
| NGC 6744 | 1.31 | 120 s × 7 | 90 s × 6 | 24 | 23.2 | 23.4 | 16 · 15 | 11 · 10 | 3 | 25.1 | 13.0 · 12.9 | V < 10.8 |

**Cómo leer la tabla:**
- **X:** masa de aire a la altura típica de observación de cada objeto.
- **S/N:** en un elemento de seeing (círculo de 4″), en el apilado, al brillo superficial medio
  del objeto. "Borde" es 2 mag/arcsec² más débil, como las partes externas de las galaxias.
- **μ límite:** brillo superficial que se detecta con S/N = 3 por elemento de seeing.
- **μ satura:** brillo superficial que satura un píxel con ese t_exp. Los núcleos más brillantes
  que ese valor saturan.
- **Estrella satura:** las estrellas del campo más brillantes que ese valor saturan en V. Es
  normal que unas pocas saturen; esas no se usan para el seeing ni para la fotometría.

**Notas por objeto:**
- **NGC 7009:** se proponen exposiciones cortas porque su interior tiene un brillo superficial
  alto. Sobra señal, así que el tiempo libre del bloque puede usarse para unas exposiciones de
  120 s que registren el halo débil.
- **NGC 1097, NGC 1291 y NGC 1316:** sus núcleos son brillantes, por eso se proponen 90 s / 60 s.
  Hay que revisar el máximo del núcleo en la prueba.
- **NGC 300 y NGC 247:** son las de más bajo brillo superficial (~23,4 mag/arcsec² en B). La S/N
  en las partes externas es baja. Si el bloque lo permite, conviene alargar las exposiciones o
  combinar ambas noches.
- **Campo:** con 0,317″/px el lado del campo es 21,6′, no 24,6′. NGC 247 (19,7′), NGC 300
  (19,4′) y la Hélice (16,3′) caben, pero dejan muy poco cielo libre para medir el fondo.

## 2. Estrellas estándar

Se usa seeing de 2,5″ y masa de aire 1,05. El t_exp propuesto es ~60 % del tiempo de saturación,
redondeado a 0,5 s.

| Estrella | B | V | t satura B | t satura V | t_exp B | t_exp V | pico B (ADU) | pico V (ADU) | S/N B | S/N V |
|---|---|---|---|---|---|---|---|---|---|---|
| HD 195500 | 7.38 | 7.32 | 4.4 s | 3.6 s | 2.5 s | 2.0 s | 22672 | 22521 | 1158 | 1154 |
| HD 202941 | 7.07 | 7.07 | 3.3 s | 2.8 s | 1.5 s | 1.5 s | 18099 | 21264 | 1033 | 1121 |
| HD 210300 | 6.59 | 6.44 | 2.1 s | 1.6 s | 1.0 s | 1.0 s | 18774 | 25325 | 1052 | 1225 |
| HD 215863 | 7.84 | 7.69 | 6.7 s | 5.0 s | 4.0 s | 2.5 s | 23747 | 20021 | 1185 | 1087 |
| HD 220881 | 7.74 | 7.45 | 6.1 s | 4.0 s | 3.5 s | 2.0 s | 22783 | 19979 | 1161 | 1086 |
| HD 562 | 7.80 | 7.66 | 6.5 s | 4.9 s | 3.5 s | 2.5 s | 21559 | 20582 | 1129 | 1103 |
| HD 8130 | 7.50 | 7.45 | 4.9 s | 4.0 s | 2.5 s | 2.0 s | 20300 | 19979 | 1095 | 1086 |
| HD 12206 | 6.81 | 6.79 | 2.6 s | 2.2 s | 1.5 s | 1.0 s | 22995 | 18346 | 1166 | 1040 |

**Recomendaciones para las estándares:**
- **Exposiciones cortas.** La cámara usa obturador electrónico (*rolling shutter*): todas las
  filas se exponen el mismo tiempo, con un desfase de 85 ms entre la primera y la última. Por
  eso 1 s es tan uniforme como 10 s. El límite lo pone el centelleo atmosférico, que agrega
  ~0,3 % de ruido en una exposición de 2 s y más en una de 1 s.
  - HD 215863, HD 562 y HD 220881 permiten 3–4 s en B.
  - HD 210300 y HD 12206 obligan a usar 1–1,5 s; con ellas conviene tomar más exposiciones.
- **Tomar 5 o más exposiciones por filtro y promediar.** La S/N de cada una ya supera 1000; la
  precisión la limitan el centelleo y el flat, no los fotones.
- **Como referencia:** en 2025 se usaron 5 s (B) y 3 s (V) para estrellas de B ≈ 8,4 y V ≈ 7,2,
  coherente con esta tabla.

## 3. Procedimiento en el telescopio (exposiciones de prueba)

En TheSkyX, al comienzo de cada bloque:

1. Apuntar, enfocar y tomar **una exposición de prueba en V** con el t_exp de la tabla.
2. Revisar el **máximo** del núcleo del objeto (o de la estrella estándar) y el **nivel de
   fondo**, ambos en ADU, y restarles el nivel de bias.
   - Si el máximo del objeto supera ~40 000 ADU: bajar t_exp.
   - Si el fondo es muy bajo (pocas decenas de ADU) y el máximo está lejos del límite: se puede
     subir t_exp, siempre que el guiado lo permita.
3. Repetir en B.
4. Anotar en la bitácora los tiempos de prueba, los valores medidos y la decisión final. El
   informe debe justificarla (punto 7 del enunciado).

Con la prueba se puede recalcular todo:
- **ZP:** el flujo F (ADU/s) de una estrella de magnitud conocida m en la prueba da
  `ZP_ADU = m + 2.5·log10(F)`. Como esta medición no está corregida por extinción, para el
  script hay que usar `ZP(e⁻/s, sobre la atmósfera) = ZP_ADU − 0,18 + k·X`.
- **Cielo:** el fondo de la prueba da la tasa de cielo, `fondo/t_exp`.
- **Nuevo cálculo:** `python tiempos_exposicion.py --zp-b … --zp-v …`.

Después de la observación, el pipeline (photutils, sin Source Extractor) mide el seeing, el
brillo del cielo y la magnitud límite reales de cada noche (`comparacion_noches.tex`). Con eso se
puede comparar lo planificado con lo obtenido.
