# Pendientes: qué confirmar, revisar y tener en cuenta

Lista de lo que falta confirmar antes y durante las noches del 8 y 15 de octubre de 2026,
y de lo que hay que considerar en el informe. Actualizada el 8/10/2026.

## Resumen

**Ya resuelto:**
- escala de placa (0,317″/px) y campo (21,6′ × 21,6′);
- datos de la cámara: ganancia, ruido de lectura, capacidad de pozo, obturador y descarga;
- hora de los FITS (UTC);
- encuadre de cada objeto en el campo;
- con los datos de 2025: zeropoint*, ruido de lectura (3,56 e⁻), nivel de bias (92,5 ADU) y seeing típico (2,5″); los tiempos de exposición ya están recalculados con ellos.

**Falta revisar en los datos de 2025:** por qué el ZP en B difiere 0,29 mag entre las dos
estándares (headers `FILTER` y `EXPTIME`).

**Falta preguntar:**  guiado, el modo de lectura, las calibraciones, el
horario y si la óptica sigue igual.

**Calibraciones y estándares:** plan con tiempos en [`calibraciones.md`](calibraciones.md). Las
magnitudes B y V de las estándares salen de Mermilliod (Johnson) con
[`estrellas_estandar.py`](estrellas_estandar.py) (astroquery).

**Falta acordar con los otros grupos:** el objeto de cada grupo, las estándares y la tabla del plan.

## 1. Antes de la observación

### Instrumento

**Resuelto** (detalles en [`verificacion_informe_2025.md`](verificacion_informe_2025.md)):

- [x] **Escala de placa y campo:** 0,317″/px y 21,6′ × 21,6′.
  - Medidos con estrellas de Gaia ([`medir_escala.py`](medir_escala.py)); coinciden con la
    solución de astrometry.net de 2025.
  - El 0,36″/px (24,6′) del informe de 2025 era incorrecto.
  - La cámara está prácticamente alineada con el norte (rotación −0,5°).
- [x] **Cámara** (ficha de Moravian y headers FITS):
  - sensor GSENSE4040 FSI, en modo de lectura 16-bit HDR;
  - ganancia de 0,85 e⁻/ADU, ruido de lectura de 3,9 e⁻ y capacidad de pozo de 56 600 e⁻;
  - el convertidor satura en 65 535 ADU, antes que el pozo;
  - obturador electrónico: las exposiciones cortas son uniformes;
  - descarga de 0,25 s por USB 3.
- [x] **Hora en los FITS:** `DATE-OBS` está en UTC (hora local = UTC − 3). El pipeline ya lo
  trata así. El informe de 2025 confundió ambas.
- [x] **Sitio (El Sauce):** −30,4597°, −70,7503°, 1600 m (header). Ya está en
  `pipeline/config.example.yaml`. Para la visibilidad, el enunciado pide Cerro Tololo.

**Medir con los datos de 2025** (notebook del pipeline: `paso_calibracion` y `paso_estandares`,
que no necesitan imágenes de ciencia):

- [x] **Corriente oscura:** no detectada. Dark − bias dio −0,039 ADU/s, es decir, los darks de
  100 s quedaron 3,9 ADU bajo el bias. Para planificar se usa 0,033 e⁻/s/px como cota
  conservadora.
- [x] **Nivel de bias.** Se suma a la señal: el margen real hasta saturar es 65 535 ADU menos el
  bias. Sale en `calibracion.json` (`bias_mediana_adu`): 92,5 ADU.
- [x] **Zeropoint:** B 21,845 y V 21,941 en ADU/s, que son 21,665 y 21,761 en e⁻/s. Queda
  0,3 mag bajo la estimación teórica: el sistema es ~25 % menos eficiente que lo supuesto.
- [ ] **Diferencia de 0,29 mag en el ZP de B entre HIP 116375 (21,99) y HIP 117678 (21,70).** En V
  coinciden a 0,01 mag. Revisar `FILTER` y `EXPTIME` en el header de una imagen B de cada
  estrella (por ejemplo, las imágenes 557 y 592). Si no hay una explicación, considerar ±0,15 mag
  de incertidumbre en el ZP de B.
- [x] **Tiempos de exposición recalculados** con esos valores, que ahora son los predeterminados
  de `tiempos_exposicion.py`. La S/N baja ~20 % y las estándares admiten ~30 % más de tiempo.
- [x] **Ruido de lectura:** 3,56 e⁻ (4,19 ADU), coherente con los 3,9 e⁻ de la ficha. Confirma el
  modo 16-bit HDR.
- [x] **Seeing de 2025:** 2,5″ de mediana (1,9–3,2″) en las 40 exposiciones de estándares.

**Preguntar al ayudante:**

- [ ] **Configuración óptica:** ¿es la misma de 2023 y 2025? La focal real es ~5860 mm; el
  `FOCALLEN = 6500` del header es solo nominal.
- [ ] **Modo de lectura en TheSkyX:** confirmar que es 16-bit HDR. En los modos de 12 bits cambian
  la ganancia y el ruido (baja ganancia: 19,5 e⁻/ADU y 34,5 e⁻).
- [ ] **Guiado:** ¿guía la montura (cámara ZWO ASI174)? Las exposiciones de 90–120 s dependen de
  eso; sin guiado, probar si 60 s salen sin estelas.
- [ ] **Temperatura del sensor:** en 2025 quedó a −14 °C sin llegar a la consigna (−16 °C); en
  2023 trabajó a −25 °C. La corriente oscura depende de la temperatura, así que los darks deben
  tomarse a la misma que la ciencia.
- [ ] **Tiempo entre exposiciones:** la descarga es de 0,25 s, pero falta el guardado en TheSkyX
  y el cambio de filtro. Se supusieron 5 s por imagen, lo que define cuántas caben en 30 min.

### Plan de observación (con los otros grupos y el ayudante)

- [ ] **Elegir el objeto de cada grupo** sin repetir.
  - Los objetos que salen tarde (NGC 1097, 1291, 1316 y 1313) van en los últimos bloques.
  - NGC 6744 va al comienzo.
  - NGC 247 y NGC 300 llenan casi todo el campo (márgenes de 1,0′ y 1,6′). Si se eligen, hay que
    asumir que el fondo se mide con poco cielo libre.
- [ ] **Elegir las estándares de cada noche** (plan en [`calibraciones.md`](calibraciones.md)):
  HD 202941 y HD 8130 al comienzo y al final, para medir el coeficiente de extinción k con dos
  masas de aire, y HD 220881 o HD 562 a mitad de la noche.
  - HD 195500 y HD 215863 se reemplazaron por HD 207480 y HD 212643, que están en Mermilliod
    (las primeras tenían B−V de Tycho, no Johnson). Avisar a los otros grupos del cambio.
  - HD 210300 y HD 12206 obligan a ~1 s y a tomar más exposiciones por el centelleo.
- [ ] **Hacer la tabla del plan** con el formato del enunciado: hora, AR, Dec, objeto, tipo,
  magnitud, distancia a la Luna y observador.
- [ ] **Horario.** Confirmar con el ayudante la ventana (21:00–01:00 o hasta las 02:00) y la
  hora de inicio (el enunciado dice desde las 19:00).
- [ ] **Calibraciones.** Plan propuesto en [`calibraciones.md`](calibraciones.md): 25 + 25 bias,
  10 darks por cada t_exp de ciencia (120, 90 y 60 s) y 10–15 flats de crepúsculo por filtro
  entre las 20:00 y las 20:35. Preguntar si los darks se pueden tomar con la cúpula cerrada y si
  hay pantalla para flats de cúpula. **Pedir darks con el mismo tiempo de exposición que la
  ciencia** (las estándares de 1–3 s se corrigen bien solo con el bias), tomados la misma noche y a la misma temperatura del sensor. En 2025 los darks de 100 s
  quedaron *bajo* el bias, así que escalar bias + corriente oscura no es confiable en este
  sensor CMOS.
- [ ] **Noche de respaldo.** Si el clima falla, se observa el jueves 22/10. Ese día la Luna está
  al ~88 %, en AR 23h20m, y sobre el horizonte toda la ventana.
  - Queda a 22° de la Hélice, 28° de NGC 247, 32° de NGC 7793, 35° de NGC 7009 y 42° de NGC 300.
  - Los objetos de AR ~3h (NGC 1097, 1291, 1316 y 1313) y NGC 6744 quedan a más de 55°.
  - El fondo será mucho más alto, sobre todo en B.

### Revisiones finales de los candidatos

- [ ] **Staralt.** Generar las curvas de visibilidad del 8/10 y del 15/10 para el objeto
  elegido y las estándares; van en el informe.
  - Staralt usa UTC−4: su marca de las 24 h corresponde a la 01:00 en Chile.
- [ ] **SIMBAD.** Revisar la ficha del objeto y de las estándares elegidas. En las nebulosas
  planetarias, SIMBAD da la magnitud de la estrella central, no la de la nebulosa: aclarar en
  el informe cuál se usa.

## 2. Durante la observación (en TheSkyX)

- [ ] **Al comienzo de la noche:** confirmar en TheSkyX el modo de lectura (16-bit HDR) y que el
  sensor llegue a la temperatura de consigna.
- [ ] **Exposición de prueba** al comienzo de cada bloque, primero en V y luego en B.
  - Revisar el **máximo** del núcleo o de la estrella y el **nivel de fondo**, ambos en ADU y
    restando el bias.
  - Si el máximo supera ~40 000 ADU, bajar el tiempo. Si el fondo es muy bajo y el máximo está
    lejos del límite, se puede subir (si el guiado lo permite).
  - **Encuadre:** revisar que el objeto quede centrado y no toque los bordes. Es crítico en
    NGC 247 (1′ libre en norte-sur) y NGC 300 (1,6′ libre en este-oeste).
- [ ] **Anotar en la bitácora:** tiempos de prueba, valores medidos, decisión final, hora (y si
  es local o UTC), filtro, número de archivo, temperatura del sensor, nubes, viento y satélites.
  El enunciado pide justificar los tiempos.
- [ ] **Enfoque.** Revisar que el FWHM de las estrellas no empeore durante la noche. En 2025 se
  midieron ~3,9″ (con la escala correcta), que es alto; puede haber influido el enfoque.
- [ ] **Satélites.** Revisar cada imagen. Con exposiciones largas, perder una cuesta más; si una
  sale contaminada, tomar otra.
- [ ] **Estándares.** Tomar 5 o más exposiciones por filtro, sin saturar. Con ~1 s tomar
  más, porque el centelleo agrega más ruido en exposiciones cortas.
- [ ] **15/10:** la Luna está alta hasta las 00:52. El fondo en B será más alto; revisarlo en la
  prueba.
- [ ] **Descargar todo** al terminar: la ciencia del grupo, las estándares y todas las
  calibraciones (bias, darks, flats y BPM).

## 3. Después, para el análisis y el informe

- [ ] **Pipeline.** Copiar `config.example.yaml` como `config.yaml` y completarlo:
  - ya trae la escala (0,317″/px), la ganancia (0,85 e⁻/ADU) y el sitio;
  - falta el objeto y sus coordenadas;
  - las estándares con sus magnitudes B y V: el bloque lo genera `estrellas_estandar.py`;
  - las rutas de los archivos de las dos noches.
- [ ] **Detección con photutils.** Se usa photutils en vez de Source Extractor. El pipeline ya
  tiene las dos configuraciones (agresiva y extendida); hay que explicar en qué se diferencian
  los parámetros de deblending.
- [ ] **Astrometría.** Usar astrometry.net (modo `archivo` o `astrometry_net`); el modo `header`
  es solo aproximado. La escala resuelta debería dar ~0,317″/px.
- [ ] **Sección "Observaciones" del informe:** detalles técnicos del telescopio y la cámara.
  - Usar la ficha de Moravian y la escala medida, no el 0,36″/px ni la focal de 6500 mm.
  - El límite de difracción es 0,22″ en B y 0,28″ en V.
- [ ] **Comparación de noches.** Usar `comparacion_noches.tex` (seeing, brillo del cielo,
  magnitud límite y Luna) para responder si un conjunto de imágenes es mejor que el otro.
- [ ] **Comparar lo planificado con lo obtenido:** S/N, cielo y seeing reales frente a los
  estimados en `tiempos_exposicion.md`.
- [ ] **Errores de 2025 que no hay que repetir** (ver `verificacion_informe_2025.md`):
  - usar `log10` y no `log`;
  - normalizar cada imagen por su propio tiempo de exposición;
  - usar magnitudes Johnson (B, V), no Tycho (BT, VT), y no cruzar las estrellas;
  - revisar que el zeropoint sea razonable (~20–23 en ADU/s, no 39);
  - calcular el error del ZP como σ/√N;
  - marcar los píxeles malos como NaN, no como 0;
  - alinear las imágenes antes de combinarlas;
  - calcular el color como B − V, cruzando las fuentes por posición;
  - poner las horas en local o en UTC, pero rotuladas correctamente;
  - usar la escala medida (0,317″/px) para el seeing y el campo.
- [ ] **Código.** Va completo y comentado en el apéndice. Todos los integrantes deben poder
  explicarlo.
