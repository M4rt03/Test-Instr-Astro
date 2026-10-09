# Pipeline Tarea 2: FIS2421 Instrumentación Astronómica

Pipeline en Python para reducir, calibrar y analizar las observaciones en B y V
del telescopio MAS500 (OMA, El Sauce). Cubre los puntos 11 a 16 del enunciado:

| Punto | Paso | Función |
|---|---|---|
| 11 | Reducción (bias, dark current, flats, alineamiento, apilado) | `paso_calibracion`, `paso_ciencia` |
| 12 | Calibración astrométrica (WCS) | `paso_astrometria` |
| 13 | Calibración fotométrica (estándares, zeropoint) | `paso_estandares` |
| 14 | Detección y catálogos con dos configuraciones + diagrama color-magnitud | `paso_catalogos` |
| 15 | Seeing por exposición, límite de difracción y variación de la PSF | `paso_seeing` |
| 16 | Imagen RGB y mapa de color B−V | `paso_color` |
| Informe | Bitácora, curvas de visibilidad, comparación entre noches | `paso_bitacora`, `paso_resumen_noche` |

Está construido a partir del material de 2025 que está en la raíz del
repositorio: se corrigieron sus errores y se agregaron las partes que faltaban
(ver más abajo).

## Instalación

```bash
cd pipeline
python -m venv .venv && source .venv/bin/activate   # opcional
pip install -r requirements.txt
```

## Uso

1. **Probar sin datos.** Esto genera un set sintético con valores verdaderos conocidos y lo procesa completo (~30 s):
   ```bash
   python run_pipeline.py --demo demo
   ```
   Los resultados quedan en `demo/resultados/`.
2. **Con tus datos.** Copia `config.example.yaml` como `config.yaml` y completa:
   - rutas (patrones glob) de bias, darks, flats, ciencia y estándares de cada noche;
   - el objeto (nombre y coordenadas) y las estándares (coordenadas y magnitudes B/V o BT/VT);
   - las imágenes malas que quieras excluir (`excluir: ["*00000529*"]`).
3. Ejecuta:
   ```bash
   python run_pipeline.py config.yaml              # todas las noches
   python run_pipeline.py config.yaml --noches noche1
   ```
   o recorre los pasos uno a uno en `notebooks/tarea2_paso_a_paso.ipynb`.
4. **Pruebas.** Comprueban que el pipeline recupera el zeropoint, el seeing, las magnitudes y los colores del set sintético:
   ```bash
   python -m pytest
   ```
5. **En el telescopio: revisión rápida de las exposiciones de prueba.** Muestra el header
   (filtro, t_exp, temperatura, ganancia y sitio), la masa de aire, el fondo del cielo, el máximo
   del centro con el t_exp máximo sin saturar, el seeing y la elipticidad. Si la imagen es de una
   estrella estándar conocida, también da su pico y el ZP:
   ```bash
   python revisar_prueba.py datos/noche1/prueba_V.fit
   python revisar_prueba.py datos/noche1/prueba.fit --estandar "HD 8130" --bias datos/noche1/Bias_1x1_00000001.fit
   python revisar_prueba.py datos/noche1/prueba_120s.fit --dark datos/noche1/Dark_120.000secs_00000001.fit
   ```
   Para medir el cielo conviene `--dark` con un dark del mismo t_exp: en 2025 el bias de 0 s quedó
   3–5 ADU sobre el nivel cero de las exposiciones.
   Conviene probarlo antes con los datos de 2025: en las estándares debería dar ZP de ~21,97–21,99
   (B, HIP 116375), ~21,70 (B, HIP 117678) y ~21,94 (V).

## Qué produce

En `resultados/<noche>/`:

| Archivo | Para el informe |
|---|---|
| `bitacora.tex`, `exposiciones.tex` | Tabla de observaciones: archivo, filtro, hora, exposición, masa de aire, seeing y desplazamiento de cada exposición |
| `figuras/visibilidad.png`, `visibilidad.tex` | Curva de visibilidad del objeto y de las estándares, con la Luna |
| `calibracion.json`, `figuras/calibracion.png` | Master bias, dark current, flats, ruido de lectura y niveles de los flats |
| `figuras/reduccion_{B,V}.png` | Imagen cruda, reducida y apilado final |
| `stack_{B,V}.fits` | Apilados en ADU/s con WCS y máscara de saturación (extensión `SATMASK`) |
| `astrometria.json` | Escala de placa, rotación y matriz CD |
| `estandares.tex`, `zeropoints.tex` | Fotometría de cada exposición de las estándares y ZP ± error por filtro |
| `seeing.json`, `figuras/seeing.png`, `psf_{B,V}.ecsv` | Seeing, comparación con la difracción y FWHM por zona del detector |
| `catalogo_{agresiva,extendida}_{B,V}.ecsv`, `catalogos.json`, `figuras/segmentacion.png` | Catálogos de las dos configuraciones y su comparación |
| `catalogo_BV.ecsv`, `figuras/color_magnitud.png` | Catálogo B y V cruzado y diagrama color-magnitud |
| `figuras/color_rgb.png`, `figuras/mapa_color_BV.png` | Imagen a color y mapa de color B−V |

En `resultados/comparacion_noches.tex` hay, para cada noche y filtro: seeing,
brillo del cielo (mag/arcsec²), magnitud límite a 5σ, zeropoint, masa de aire e
iluminación y distancia de la Luna. Es la base para la discusión
"¿hay un conjunto de imágenes mucho mejor que el otro?".

## Qué se corrigió respecto del material de 2025

| Problema en el código anterior | Cómo queda aquí |
|---|---|
| Las magnitudes del catálogo usaban `np.log` (logaritmo natural) mientras que el ZP se calculó con `log10`; salían magnitudes fuera de escala y B−V entre −10 y +10 | `log10` en todas partes; hay una prueba que verifica B−V contra valores verdaderos |
| El índice de color se calculaba como `V − B` y se emparejaban estrellas por su orden en la tabla | `B − V`, con cruce por posición y emparejamiento uno a uno (`detection.crossmatch`) |
| Todas las estándares se dividían por el `EXPTIME` de una sola imagen (5 s), aunque las de V eran de 3 s | Cada imagen se normaliza por su propio tiempo (todo queda en ADU/s) |
| Etiquetas de HIP 116375 y HIP 117678 cruzadas; magnitudes Tycho usadas como si fueran Johnson | Magnitudes por estrella en el YAML y conversión Tycho → Johnson (ESA 1997) |
| Error del ZP = σ/N | σ/√N, con una medición por exposición (no una sola) |
| Seeing medido en una sola estrella estándar (posiblemente saturada); la tabla tenía errores de cálculo | Ajuste gaussiano 2D a decenas de estrellas aisladas y no saturadas de cada exposición |
| Imágenes combinadas sin alinear | Alineamiento sub-píxel por correlación cruzada; B alineado a V |
| BPM multiplicada (píxeles malos = 0) | Píxeles malos = NaN, que se ignoran en las estadísticas |
| Flats combinados sin normalizar cada uno | Cada flat se normaliza por su propia mediana antes de combinar |
| Un solo `detect_threshold` global en una imagen con nebulosa | Fondo 2D, filtro gaussiano y umbral en sigmas de la imagen filtrada |
| Una sola configuración de detección (se pedían dos) | `agresiva` y `extendida`, con parámetros documentados |
| Flujo = suma del segmento (incluye la nebulosa y las vecinas) | Apertura de 1 FWHM sobre la imagen sin nebulosa, fondo local y corrección de apertura con una PSF empírica |
| Nombres de archivo escritos a mano (decenas de líneas repetidas) | Patrones glob en un YAML y bucles |
| Sin comparación entre noches ni límite de difracción | Tabla comparativa por noche y cociente seeing / difracción |

## Lo que tienes que revisar con tus datos

- **Ganancia y saturación** (`instrumento.ganancia_e_adu`, `saturacion_adu`): revisa el header (`EGAIN`, `GAIN`) o el manual de la cámara. La ganancia solo afecta a los errores.
- **Coordenadas del sitio**: verifícalas con el manual del MAS500. Solo se usan si el header no trae `AIRMASS` y para las curvas de visibilidad.
- **`DATE-OBS`** se asume en UTC, que es el estándar FITS.
- **Extinción**: los valores `k` por defecto son típicos de Cerro Tololo. Si observas estándares a distintas masas de aire, puedes ajustar `k` con las tablas de `estandares.ecsv`.
- **Astrometría**: el modo `header` solo da un WCS aproximado. Para el informe usa astrometry.net (modo `archivo` o `astrometry_net`).
- **Imágenes con satélites**: revisa las exposiciones individuales (`exposiciones.ecsv` muestra el fondo y el seeing de cada una). El sigma-clipping elimina la mayoría de las trazas, pero puedes excluir una imagen con `excluir`.
- **Estrellas débiles sobre la parte brillante de la nebulosa**: tienen errores sistemáticos mayores (el fondo tiene estructura a la escala de la PSF). El error reportado incluye un término por esa estructura, pero conviene mencionarlo en la discusión.
- **Nombres de filtros**: el pipeline compara el `FILTER` del header con `B`/`V` y avisa si no coinciden. Si tu cámara escribe otro nombre (por ejemplo `Johnson B`), ajusta las claves del YAML.
