# Candidatos para la Tarea 2 (2026): objetos y estrellas estándar

Listas para elegir el objeto de cada grupo y las estrellas estándar comunes
para las noches de observación del **jueves 8 y el jueves 15 de octubre de 2026**
con el telescopio MAS500 (campo de 24,6′ × 24,6′, filtros B y V).

## Criterios del enunciado

- El objeto debe ser una **nebulosa o galaxia**, con magnitud aparente entre **~6 y 10**
  en la banda correspondiente.
- Debe ser **visible ambas noches antes de la 01:00** desde Cerro Tololo.
- Las estándares deben tener magnitud **entre 6 y 8 en cada banda** (B y V).
- Criterio adicional: **distancia a la Luna la noche del 15/10** (creciente al 26 %).

## Cómo se calculó

- **Sitio:** Cerro Tololo (−30,169°, −70,807°, 2200 m). La hora local es UTC−3.
- **Ventana útil:** desde el fin del crepúsculo náutico (Sol a −12°) hasta la 01:00.
  - 8/10: puesta del Sol 19:49, fin del crepúsculo náutico 20:42 y astronómico 21:11.
  - 15/10: puesta del Sol 19:54, fin del crepúsculo náutico 20:47 y astronómico 21:17.
- **"h > 30°"** son las horas de esa ventana con altura mayor a 30° (masa de aire < 2).
- **"Alt. máx"** es la altura máxima dentro de la ventana, no la de culminación.
- **Fuentes:**
  - Objetos: [OpenNGC](https://github.com/mattiaverga/OpenNGC) y RC3 (VizieR VII/155).
  - Estrellas: Hipparcos (VizieR I/239), con magnitudes y tipos espectrales de SIMBAD.
- **Verificación de las estrellas:** son estrellas simples en SIMBAD, sin bandera de
  variabilidad (`VarFlag`) ni solución de componentes (`MultFlag`) en Hipparcos.
- **Verificación de las curvas (5/10/2026):** se recalcularon en forma independiente las
  alturas, minuto a minuto, con las fórmulas de Meeus (*Astronomical Algorithms*) para el
  Sol y la Luna, incluyendo la paralaje lunar. Coinciden con el cálculo original con
  astropy dentro de 0,1 h y 2°. Para la figura del informe, usar Staralt (ver al final)
  o `paso_bitacora` del pipeline. Staralt (catserver.ing.iac.es) y JPL Horizons no fueron
  accesibles desde el entorno donde se hizo la verificación.
- **Verificación en SIMBAD y RC3 (5/10/2026):**
  - Las 8 estándares coinciden con SIMBAD en coordenadas (< 1″), B, V y tipo espectral, y
    SIMBAD las clasifica como estrellas simples (`*`), no como variables ni dobles.
  - Las coordenadas de la versión anterior diferían de SIMBAD en hasta 12″ (los segundos
    se habían truncado). Ahora todas las coordenadas son las de SIMBAD.
  - Las magnitudes se corrigieron donde hacía falta (ver "Origen de las magnitudes").

## Condiciones de las noches

| Noche | Luna | Consecuencia |
|---|---|---|
| Jue 8/10/2026 | Casi nueva (3 % iluminada), bajo el horizonte toda la noche | Cielo oscuro toda la ventana |
| Jue 15/10/2026 | Creciente (26 %), en AR 17h31m, Dec −28° (límite Ofiuco–Sagitario). Está a 48° de altura a las 20:47 y se pone a las 00:51 | Los objetos de AR 18h–20h quedan cerca de la Luna. Después de las 00:51 el cielo queda oscuro |

**Distancias a la Luna el 15/10 a las 20:47** (la Luna avanza ~0,5°/h hacia el este):

- **< 20°, descartados:** M8 y M20 (9°) y M17 (17°).
- **30°–40°, aceptables:**
  - estándares HIP 97210 (30°) y HIP 98926 (35°);
  - objetos NGC 6818 y NGC 6822 (34°) y NGC 6744 (39°).

  La Luna está al 26 %, así que el fondo de cielo sube pero sigue siendo usable. Si es
  posible, conviene observarlos al final de su ventana o elegir otra opción.
- **> 50°, sin problema:** todo lo demás.

## 1. Objetos candidatos (10)

| # | Objeto | Tipo | AR (J2000) | Dec (J2000) | B | V | Tamaño | Alt. máx | h > 30° 8/10 · 15/10 | Dist. Luna 15/10 | Notas |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **NGC 253** (galaxia del Escultor) | galaxia espiral | 00:47:33.1 | −25:17:19.7 | 8,0 | 7,2 | 27′ × 5′ | 85° | 4,2 · 4,2 | 94° | Brillante y con mucha estructura (polvo, brotes de formación estelar). Estándar cercana: HIP 4442 |
| 2 | **NGC 7293** (nebulosa de la Hélice) | nebulosa planetaria | 22:29:38.5 | −20:50:13.7 | 7,5 | 7,3 | 16′ | 81° | 4,3 · 4,2 | 68° | Cabe completa en el campo. Brillo superficial bajo: requiere exposiciones largas |
| 3 | **NGC 55** | galaxia de canto | 00:14:53.6 | −39:11:47.9 | 8,6 | 7,9 | 30′ × 3′ | 81° | 4,3 · 4,2 | 81° | Más larga que el lado del campo: cabe orientada en la diagonal |
| 4 | **NGC 300** | galaxia espiral | 00:54:53.4 | −37:41:03.2 | 8,7 | 8,1 | 19′ × 13′ | 83° | 4,3 · 4,2 | 88° | Brillo superficial bajo |
| 5 | **NGC 7793** | galaxia espiral | 23:57:49.8 | −32:35:27.7 | 9,6 | 9,1 | 10′ × 6′ | 88° | 4,3 · 4,2 | 81° | Pasa casi por el cenit; tamaño cómodo para el campo |
| 6 | **NGC 7009** (nebulosa de Saturno) | nebulosa planetaria | 21:04:10.8 | −11:21:48.6 | 8,3 | 8,0 | 0,7′ | 71° | 4,3 · 4,2 | 53° | Pequeña y de alto brillo superficial; el campo tiene muchas estrellas para medir el seeing |
| 7 | **NGC 246** (nebulosa de la Calavera) | nebulosa planetaria | 00:47:03.3 | −11:52:19.0 | 8,0 | 10,9 | 4′ | 72° | 3,8 · 4,2 | 101° | V en el límite del rango; comparte estándar con NGC 253 |
| 8 | **NGC 247** | galaxia espiral | 00:47:08.6 | −20:45:37.4 | 9,7 | 9,1 | 21′ × 7′ | 80° | 4,1 · 4,2 | 96° | A 4,5° de NGC 253 y a 5° de HIP 4442. Brillo superficial bajo (SIMBAD la clasifica como LSB). Reemplaza a M17 (ver abajo) |
| 9 | **NGC 6744** | galaxia espiral | 19:09:46.1 | −63:51:26.9 | 9,1 | 8,3 | 16′ × 10′ | 55° | 4,3 · 4,0 | 39° | Circumpolar; mejor al comienzo de la noche. El 15/10 la Luna está a 39° |
| 10 | **NGC 6818** (Little Gem) | nebulosa planetaria | 19:43:57.8 | −14:09:13.4 | 9,9 | 9,3 | 0,8′ | 73° | 3,8 · 3,3 | 34° | Observar temprano (baja de 30° a las 00:30 el 8/10 y a las 00:02 el 15/10). Su alto brillo superficial tolera la Luna |

**Reemplazos posibles:**
- NGC 1097 (galaxia; B 10,2, V 9,5; 9′ × 6′). Sale tarde: supera 30° desde las 22:38 (8/10) y las
  22:10 (15/10), y está a 110° de la Luna.
- NGC 6822 (galaxia de Barnard; B 9,3; 15′). No hay una V confiable en los catálogos consultados. Tiene brillo superficial muy bajo y
  está a 34° de la Luna el 15/10: solo si se acepta un fondo más alto esa noche.

**No recomendados por la Luna del 15/10:**
- **NGC 6618 (M17, objeto del informe de 2025):** queda a 17° de la Luna. Además solo está
  sobre 30° hasta las 23:11 (8/10) y las 22:44 (15/10), es decir 2,5 y 2,0 h.
- **M8 (Laguna) y M20 (Trífida):** quedan a 9° de la Luna, y M8 es más grande que el campo.

### Origen de las magnitudes

- **NGC 253 y NGC 7793:** RC3, con V = B_T − (B−V)_T. OpenNGC da V = 11,1 para NGC 253, que es un error.
- **NGC 247 y NGC 1097:** RC3. SIMBAD da casi lo mismo: B 9,61 y V 9,10 para NGC 247; B 9,97 y V 9,48 para NGC 1097.
- **NGC 55, NGC 300 y la V de NGC 6744 y NGC 6818:** SIMBAD.
- **B de NGC 6744:** RC3 da 9,14 y SIMBAD 9,28.
- **NGC 7793:** SIMBAD da B 9,74 y V 9,28, unas 0,15 mag más débil que RC3.
- **NGC 6822:** B de RC3. SIMBAD da V 8,1 junto con B 18,0, un dato claramente erróneo, así que su V no es confiable.
- **Nebulosas planetarias:** magnitudes **integradas** de OpenNGC. SIMBAD les asigna la
  magnitud de la estrella central (B ≈ 11,5–13,5), que **no** es la de la nebulosa. Hay que
  aclarar en el informe cuál se usa.

## 2. Estrellas estándar candidatas (8)

Son enanas de tipo A0–A9 cerca del cenit (δ ≈ −22° a −31°), repartidas en ascensión recta
para que siempre haya una alta a cualquier hora de la ventana.

| # | Estrella | HD | AR (J2000) | Dec (J2000) | V | B | B−V | Tipo | Alt. máx | h > 30° 8/10 · 15/10 | Dist. Luna 15/10 | Sirve para |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **HIP 97210** | 186417 | 19:45:24.3 | −30:54:10.8 | 6,81 | 6,93 | 0,12 | A1V | 85° | 4,3 · 3,8 | 30° | Comienzo de noche el 8/10: NGC 6744, NGC 6818 |
| 2 | **HIP 98926** | 190285 | 20:05:12.0 | −26:48:47.8 | 7,22 | 7,31 | 0,09 | A1V | 87° | 4,3 · 4,0 | 35° | Comienzo de noche el 8/10 |
| 3 | **HIP 105318** | 202941 | 21:19:51.4 | −27:12:32.9 | 7,05 | 7,05 | 0,00 | A0V | 87° | 4,3 · 4,2 | 51° | Estándar temprana del 15/10; NGC 7009 |
| 4 | **HIP 107766** | 207480 | 21:49:53.9 | −27:24:14.2 | 7,11 | 7,17 | 0,06 | A1V | 87° | 4,3 · 4,2 | 57° | Estándar temprana del 15/10; NGC 7009, NGC 7293 |
| 5 | **HIP 110861** | 212852 | 22:27:39.3 | −26:25:07.7 | 7,20 | 7,41 | 0,21 | A9V | 86° | 4,3 · 4,2 | 65° | NGC 7293 (a 6°) |
| 6 | **HIP 116750** | 222332 | 23:39:41.6 | −22:32:00.2 | 7,29 | 7,41 | 0,12 | A0V | 82° | 4,3 · 4,2 | 82° | NGC 7793, NGC 55 |
| 7 | **HIP 4442** | 5524 | 00:56:49.1 | −25:21:48.1 | 7,24 | 7,33 | 0,09 | A2/3V | 85° | 4,1 · 4,2 | 96° | NGC 253 (a 2°), NGC 247, NGC 246, NGC 300 |
| 8 | **HIP 6897** | 9063 | 01:28:47.9 | −24:47:52.5 | 7,04 | 7,27 | 0,23 | A8V | 81° | 3,5 · 4,0 | 102° | Final de la noche |

Las magnitudes B y V son las de SIMBAD. Las 8 cumplen 6 ≤ B, V ≤ 8.

**Descartadas en la verificación:**
- **HIP 112362 y HIP 6393:** Hipparcos las marca como posibles variables (`VarFlag = 1`).
- **HIP 4496 y HIP 6507:** Hipparcos las resuelve como dobles (`MultFlag = C`).
- **HIP 116375 y HIP 117678** (las del informe de 2025): son gigantes K con B ≈ 8,1, fuera del
  rango 6–8 en B.

## Recomendaciones para el plan de observación

- **Estándar temprana y estándar tardía:** así se observan a masas de aire distintas y se
  puede **medir** el coeficiente de extinción k. Si no, el pipeline usa el valor típico de
  Cerro Tololo (`extincion` en `config.yaml`).
  - Como el enunciado pide el mismo plan ambas noches, la opción más simple es
    **HIP 105318 o HIP 107766** al comienzo (a ≥ 51° de la Luna el 15/10) y
    **HIP 4442 o HIP 6897** al final.
  - HIP 97210 y HIP 98926 sirven el 8/10, pero el 15/10 quedan a 30°–35° de la Luna.
    Son estrellas brillantes y con exposiciones cortas el efecto es pequeño, pero conviene
    evitarlas.
- **Tiempos de exposición de las estándares:** con V ≈ 7, conviene partir las pruebas en 1–3 s
  y revisar que el pico de la estrella quede bajo ~60 000 ADU. El pipeline excluye del
  zeropoint las mediciones saturadas.
- **Verificación en SIMBAD:** antes de fijar el plan, conviene revisar la ficha de cada
  estrella elegida (https://simbad.cds.unistra.fr) por si hay información nueva.

## Curvas de visibilidad con Staralt

En [Staralt](http://catserver.ing.iac.es/staralt/):

1. Elegir el modo **Staralt**, la fecha **8 de octubre de 2026** (después **15 de octubre**) y
   el observatorio **Cerro Tololo Observatory (Chile)**. Staralt usa la fecha del inicio de la noche.
2. Pegar la lista de coordenadas siguiente (formato `nombre hh mm ss ±dd mm ss`).
3. Dejar marcada la opción que muestra la Luna. La leyenda indica la distancia a la Luna de cada objeto.

```
NGC253   00 47 33 -25 17 20
NGC7293  22 29 38 -20 50 14
NGC55    00 14 54 -39 11 48
NGC300   00 54 53 -37 41 03
NGC7793  23 57 50 -32 35 28
NGC7009  21 04 11 -11 21 49
NGC246   00 47 03 -11 52 19
NGC247   00 47 09 -20 45 37
NGC6744  19 09 46 -63 51 27
NGC6818  19 43 58 -14 09 13
HIP97210  19 45 24 -30 54 11
HIP98926  20 05 12 -26 48 48
HIP105318 21 19 51 -27 12 33
HIP107766 21 49 54 -27 24 14
HIP110861 22 27 39 -26 25 08
HIP116750 23 39 42 -22 32 00
HIP4442   00 56 49 -25 21 48
HIP6897   01 28 48 -24 47 52
```

Valores esperados para comparar con las figuras de Staralt (altura en grados):

| Objeto | 20:30 | 21:30 | 22:30 | 23:30 | 00:30 | 01:00 |
|---|---|---|---|---|---|---|
| NGC 253 (8/10 · 15/10) | 26 · 32 | 39 · 45 | 52 · 58 | 65 · 71 | 78 · 83 | 83 · 85 |
| NGC 7293 | 54 · 60 | 67 · 72 | 78 · 80 | 79 · 75 | 68 · 63 | 62 · 56 |
| NGC 55 | 37 · 42 | 49 · 54 | 60 · 66 | 72 · 76 | 80 · 81 | 81 · 78 |
| NGC 300 | 29 · 34 | 41 · 46 | 53 · 58 | 64 · 70 | 76 · 80 | 81 · 83 |
| NGC 7793 | 39 · 45 | 51 · 57 | 64 · 70 | 77 · 82 | 88 · 84 | 83 · 78 |
| NGC 7009 | 65 · 69 | 71 · 71 | 68 · 64 | 58 · 53 | 46 · 40 | 40 · 34 |
| NGC 246 | 21 · 27 | 34 · 40 | 46 · 52 | 58 · 63 | 68 · 71 | 71 · 71 |
| NGC 247 | 25 · 30 | 37 · 43 | 50 · 56 | 63 · 69 | 75 · 79 | 79 · 80 |
| NGC 6744 | 56 · 54 | 52 · 50 | 48 · 45 | 41 · 39 | 35 · 32 | 32 · 29 |
| NGC 6818 | 74 · 72 | 67 · 62 | 56 · 50 | 43 · 37 | 30 · 24 | 24 · 18 |
| Luna (solo 15/10) | 51 | 39 | 27 | 15 | 4 | −2 (se pone 00:51) |
