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
medio, `μ = m + 2.5 log10(área en arcsec²)`, y la S/N se calcula en un círculo de 4″ de
diámetro (~1,6 veces el seeing típico).

### Saturación

- **Estrella:** el píxel central recibe una fracción `1/(1.133·FWHM_px²)` del flujo total.
  Satura cuando esa fracción, sumada al cielo, supera el límite de linealidad.
- **Objeto extendido:** un píxel satura cuando `10^(0.4·(ZP − μ))·escala²·t` llega a ese límite.

### Parámetros usados

| Parámetro | Valor | Origen |
|---|---|---|
| Escala de placa | 0,317″/px | Medida con estrellas de Gaia ([`medir_escala.py`](medir_escala.py): 0,3168″/px); coincide con astrometry.net en 2025 (0,3163″/px). El 0,36″/px que usó el informe de 2025 era incorrecto |
| Ganancia | 0,85 e⁻/ADU | Header `GAIN` de los datos de 2025 |
| Ruido de lectura | 3,56 e⁻ | **Medido** en los bias de 2025 (4,19 ADU). La ficha de Moravian da 3,9 e⁻ para el modo 16-bit HDR |
| Corriente oscura | 0,033 e⁻/s/px | Cota conservadora. En 2025 **no se detectó**: dark − bias dio −0,039 ADU/s (ver debajo). Moravian no la publica |
| Límite de saturación | 40 000 ADU sobre el bias (medido: 92,5 ADU) | Conservador. El ADC satura en 65 535 ADU, antes que el pozo (56 600 e⁻ = 66 600 ADU); según Moravian, el modo HDR es lineal en todo el rango |
| ZP sobre la atmósfera | B 21,665 · V 21,761 (e⁻/s) | **Medido** con las estándares de 2025 (`zeropoints.ecsv` − 0,18); ver detalle debajo |
| Extinción | k_B = 0,25 · k_V = 0,15 | Valores típicos de Cerro Tololo |
| Cielo sin Luna (8/10) | B 22,7 · V 21,8 mag/arcsec² | Valores típicos de Cerro Tololo |
| Cielo con Luna al 25 % (15/10) | B 22,1 · V 21,6 mag/arcsec² | Valores típicos para esa fase |
| Seeing | 2,5″ típico · 2,0″ para la saturación | **Medido** por el pipeline en las 40 exposiciones de estándares de 2025: mediana 2,5″, rango 1,9–3,2″. El informe de 2025 daba 3,9–4,4″ por usar mal la escala y el método |

**El ZP medido frente al estimado:**
- La estimación teórica daba B 21,98 y V 22,05. Se calculó con un área efectiva de 1570 cm²
  (espejo de 0,5 m con una obstrucción de ~45 %), el flujo de Vega en B y V (Bessell et al. 1998)
  y una eficiencia total de 0,30 en B y 0,48 en V.
- El ZP medido con HIP 116375 y HIP 117678 (2025) es 0,32 mag más bajo en B y 0,29 en V. Es
  decir, el sistema real (espejos, filtros y eficiencia cuántica) es un ~25 % menos eficiente que
  lo supuesto, dentro de la incertidumbre de ±0,5 mag de la estimación.
- Las tablas usan los valores medidos.
- **En V las dos estrellas coinciden** (21,942 y 21,941, con 0,01 mag de dispersión).
- **En B difieren en 0,29 mag:** HIP 116375 da 21,989 y HIP 117678 da 21,701, cada una con solo
  0,008 mag de dispersión entre sus 10 exposiciones.
  - Según los catálogos, HIP 117678 debería verse un 16 % más brillante que HIP 116375 en B, pero
    se midió un 11 % más débil.
  - Las magnitudes de catálogo coinciden entre Hipparcos, Tycho-2 y SIMBAD, y la diferencia de
    masa de aire (0,13) o de color (0,1 mag) no alcanza a explicarla.
  - Lo más probable es un problema de las imágenes en B de una de las dos estrellas: otro filtro
    en la rueda, otro tiempo de exposición real o nubes durante esa serie. Hay que revisar el
    header (`FILTER`, `EXPTIME`) de una imagen B de cada estrella.
  - Para planificar se usa el promedio (21,845 en ADU/s). Un error de ±0,15 mag cambia la S/N
    solo en ±12 %. Para la fotometría del informe, en cambio, sería un error sistemático
    importante en B y en B−V.

**Corriente oscura negativa:** los darks de 100 s quedaron 3,9 ADU *por debajo* del bias. No es
físico. Indica que la corriente oscura es menor de lo que se puede medir así y que el nivel de
bias no es idéntico entre una exposición de 0 s y una de 100 s (algo común en sensores CMOS).
Consecuencia práctica: en las noches de 2026 conviene tomar **darks con el mismo tiempo de
exposición que la ciencia** (y que las estándares) y restarlos directamente, en lugar de usar
bias + corriente oscura escalada.

### Cómo medir los valores reales con los datos de 2025

Con el notebook `pipeline/notebooks/tarea2_paso_a_paso.ipynb`, ejecutando solo
`paso_calibracion` y `paso_estandares` (no necesitan las imágenes de ciencia), se obtiene:

- **`calibracion.json`:** `ruido_lectura_e`, `dark_current_mediana_adu_s` (multiplicar por la
  ganancia para pasar a e⁻/s) y `bias_mediana_adu`.
- **`zeropoints.ecsv`:** el ZP en ADU/s, ya corregido por extinción. Para el script, el ZP en
  e⁻/s es `ZP_ADU + 2.5·log10(0,85)` ≈ `ZP_ADU − 0,18`.
- **`estandares.ecsv`:** el FWHM y el ZP de cada exposición, para `--seeing` y para revisar
  que las estrellas no estén saturadas.

Luego: `python tiempos_exposicion.py --zp-b … --zp-v … --ruido … --oscuridad … --seeing …`

## Resultado principal: dominan el ruido de lectura y la corriente oscura

Con píxeles de 0,317″, cada píxel recibe muy poca luz del cielo: entre 0,03 y 0,10 e⁻/s.
Para que el ruido del cielo domine sobre el ruido de lectura (cielo > 10·RN² por píxel) harían
falta exposiciones de 20 a 70 minutos. Por lo tanto:

- **Conviene tomar pocas exposiciones largas en lugar de muchas cortas.** Con el mismo
  tiempo total (~14 min en B, NGC 300):

  | Exposición | S/N |
  |---|---|
  | 30 s × 28 | 7,0 |
  | 60 s × 14 | 9,2 |
  | 120 s × 7 | 11,6 |
  | 300 s × 3 | 14,8 |

- **El límite práctico de t_exp lo ponen otras cosas:** el seguimiento o guiado de la montura,
  los satélites (con pocas exposiciones, perder una cuesta más) y la saturación del núcleo del
  objeto. Por eso se propone 90–120 s y no 300 s. Si el guiado funciona bien, se puede subir.
- **La Luna del 15/10 casi no cambia la S/N** (NGC 300 en B: 11,6 → 10,9), porque el cielo no
  es la fuente dominante de ruido. Sí cambia el nivel de fondo, que debe sustraerse bien.
- **Los 30 s de 2025 (M17) quedaban en el régimen dominado por el ruido de lectura.** Para M17 no
  importaba, porque es muy brillante. Para galaxias de ~23 mag/arcsec² sí importa.

## 1. Objetos: plan por bloque de 30 min

Los tiempos suman ~24 min de exposición más lecturas, lo que deja ~6 min para apuntar, enfocar
y tomar las pruebas. B recibe más tiempo porque el cielo es más oscuro y la eficiencia, menor.

| Objeto | X | B: t_exp × N | V: t_exp × N | Tiempo total (min) | μB medio | μV medio | S/N B 8/10 · 15/10 | S/N V 8/10 · 15/10 | S/N B borde (μ+2) 8/10 | μ límite B 8/10 (S/N=3) | μ satura B · V | Estrella satura V |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| NGC 7293 | 1.04 | 120 s × 7 | 90 s × 6 | 24 | 22.2 | 22.0 | 34 · 32 | 32 · 31 | 6 | 24.9 | 12.8 · 12.7 | V < 11.0 |
| NGC 7009 | 1.06 | 30 s × 12 | 20 s × 12 | 12 | 16.2 | 15.9 | 739 · 738 | 761 · 760 | 269 | 23.9 | 11.3 · 11.0 | V < 9.4 |
| NGC 300 | 1.06 | 120 s × 7 | 90 s × 6 | 24 | 23.4 | 23.3 | 12 · 11 | 10 · 10 | 2 | 24.9 | 12.8 · 12.7 | V < 11.0 |
| NGC 7793 | 1.04 | 120 s × 7 | 90 s × 6 | 24 | 22.8 | 22.4 | 20 · 19 | 22 · 22 | 3 | 24.9 | 12.8 · 12.7 | V < 11.0 |
| NGC 247 | 1.06 | 120 s × 7 | 90 s × 6 | 24 | 23.4 | 22.9 | 12 · 11 | 15 · 14 | 2 | 24.9 | 12.8 · 12.7 | V < 11.0 |
| NGC 1097 | 1.15 | 90 s × 9 | 60 s × 9 | 24 | 23.3 | 23.0 | 12 · 11 | 12 · 12 | 2 | 24.8 | 12.4 · 12.2 | V < 10.6 |
| NGC 1291 | 1.22 | 90 s × 9 | 60 s × 9 | 24 | 23.1 | 22.4 | 13 · 13 | 19 · 19 | 2 | 24.8 | 12.4 · 12.2 | V < 10.6 |
| NGC 1316 | 1.22 | 90 s × 9 | 60 s × 9 | 24 | 23.1 | 22.2 | 14 · 13 | 24 · 24 | 2 | 24.8 | 12.4 · 12.2 | V < 10.6 |
| NGC 1313 | 1.41 | 120 s × 7 | 90 s × 6 | 24 | 23.3 | 23.1 | 12 · 11 | 11 · 11 | 2 | 24.9 | 12.7 · 12.6 | V < 11.0 |
| NGC 6744 | 1.31 | 120 s × 7 | 90 s × 6 | 24 | 23.2 | 23.4 | 14 · 13 | 9 · 9 | 2 | 24.9 | 12.7 · 12.6 | V < 11.0 |

**Cómo leer la tabla:**
- **X:** masa de aire a la altura típica de observación de cada objeto.
- **S/N:** en un círculo de 4″ de diámetro, en el apilado, al brillo superficial medio
  del objeto. "Borde" es 2 mag/arcsec² más débil, como las partes externas de las galaxias.
- **μ límite:** brillo superficial que se detecta con S/N = 3 en ese círculo.
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
  en las partes externas es baja (S/N ~11 al brillo medio en B). Si hay guiado, conviene subir a
  180 s × 5 en B (S/N 13,3 con el mismo tiempo total) o combinar ambas noches.
- **Campo:** con 0,317″/px el lado del campo es 21,6′, no 24,6′. NGC 247 (19,7′), NGC 300
  (19,4′) y la Hélice (16,3′) caben, pero dejan muy poco cielo libre para medir el fondo.

## 2. Estrellas estándar

Se usa seeing de 2,0″ (el mejor observado en 2025) y masa de aire 1,05. El t_exp propuesto es ~60 % del tiempo de saturación,
redondeado a 0,5 s.

| Estrella | B | V | t satura B | t satura V | t_exp B | t_exp V | pico B (ADU) | pico V (ADU) | S/N B | S/N V |
|---|---|---|---|---|---|---|---|---|---|---|
| HD 195500 | 7.38 | 7.32 | 3.8 s | 3.0 s | 2.0 s | 1.5 s | 21203 | 20224 | 899 | 877 |
| HD 202941 | 7.07 | 7.07 | 2.8 s | 2.4 s | 1.5 s | 1.0 s | 21157 | 16973 | 898 | 803 |
| HD 210300 | 6.59 | 6.44 | 1.8 s | 1.3 s | 1.0 s | 1.0 s | 21947 | 30323 | 914 | 1076 |
| HD 215863 | 7.84 | 7.69 | 5.8 s | 4.2 s | 3.0 s | 2.5 s | 20821 | 23972 | 890 | 956 |
| HD 220881 | 7.74 | 7.45 | 5.3 s | 3.3 s | 3.0 s | 2.0 s | 22829 | 23922 | 933 | 955 |
| HD 562 | 7.80 | 7.66 | 5.6 s | 4.1 s | 3.0 s | 2.0 s | 21602 | 19715 | 907 | 866 |
| HD 8130 | 7.50 | 7.45 | 4.2 s | 3.3 s | 2.5 s | 2.0 s | 23731 | 23922 | 951 | 955 |
| HD 12206 | 6.81 | 6.79 | 2.2 s | 1.8 s | 1.0 s | 1.0 s | 17921 | 21967 | 826 | 915 |

**Recomendaciones para las estándares:**
- **Exposiciones cortas.** La cámara usa obturador electrónico (*rolling shutter*): todas las
  filas se exponen el mismo tiempo, con un desfase de 85 ms entre la primera y la última. Por
  eso 1 s es tan uniforme como 10 s. El límite lo pone el centelleo atmosférico, que agrega
  ~0,3 % de ruido en una exposición de 2 s y más en una de 1 s.
  - HD 215863, HD 562 y HD 220881 permiten ~3 s en B y 2–2,5 s en V.
  - HD 210300 y HD 12206 obligan a usar ~1 s; con ellas conviene tomar más exposiciones.
- **Tomar 5 o más exposiciones por filtro y promediar.** La S/N de cada una ya es de ~800–1000; la
  precisión la limitan el centelleo y el flat, no los fotones.
- **Como referencia:** en 2025 se usaron 5 s (B) y 3 s (V) para estrellas de B ≈ 8,1 y V ≈ 7,1,
  con seeing de ~2,5″, y ninguna exposición saturó. Es coherente con esta tabla: estas estándares
  son ~0,5–1 mag más brillantes en B, por eso sus tiempos son más cortos.

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
