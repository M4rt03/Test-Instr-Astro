# Pendientes: qué confirmar, revisar y tener en cuenta

Lista de lo que falta confirmar antes y durante las noches del 8 y 15 de octubre de 2026,
y de lo que hay que considerar en el informe.

## 1. Antes de la observación

### Instrumento (preguntar al ayudante o medir con los datos de 2025)

- [x] **Escala de placa y campo.** Resuelto: 0,317″/px y un campo de 21,6′ × 21,6′, medidos con
  estrellas de Gaia ([`medir_escala.py`](medir_escala.py)) y coincidentes con astrometry.net en 2025.
  Ver [`verificacion_informe_2025.md`](verificacion_informe_2025.md).
  - NGC 247 (19,7′), NGC 300 (19,4′) y la Hélice (16,3′) caben, pero dejan muy poco cielo libre
    para el fondo. NGC 253 (26,8′) y NGC 55 (29,9′) no caben.
  - Confirmar con el ayudante que la configuración óptica sigue siendo la misma de 2023 y 2025.
- [ ] **Ruido de lectura y corriente oscura.** Se supusieron 3,7 e⁻ y 0,05 e⁻/s/px. Medirlos
  corriendo `paso_calibracion` en el notebook del pipeline con los bias y darks de 2025
  (`calibracion.json`).
- [ ] **Zeropoint.** Se estimó en teoría (±0,5 mag). Medirlo con `paso_estandares` y las
  estándares de 2025 (`zeropoints.ecsv`). Luego recalcular:
  `python tiempos_exposicion.py --zp-b … --zp-v … --ruido … --oscuridad …`.
- [ ] **Ganancia.** Confirmar que sigue siendo 0,85 e⁻/ADU (header `GAIN` en 2023 y 2025).
- [ ] **Temperatura del sensor.** En 2025 el sensor quedó a −14 °C sin llegar a la consigna
  (−16 °C); en 2023 trabajó a −25 °C. La corriente oscura depende de la temperatura: tomar los
  darks a la misma temperatura que la ciencia.
- [ ] **Hora en los FITS.** `DATE-OBS` está en UTC (hora local = UTC − 3). El informe de 2025
  confundió ambas.
- [ ] **Nivel de bias y saturación.** El nivel del bias se suma a la señal: el límite útil es
  ~65 535 ADU menos el bias. Se usó un límite conservador de 40 000 ADU sobre el bias.
- [ ] **Guiado.** Preguntar si la montura guía (cámara ZWO ASI174). Las exposiciones de
  90–120 s dependen de eso; sin guiado, probar si 60 s salen sin estelas.
- [ ] **Tiempo de descarga y de cambio de filtro.** Se supusieron 5 s por imagen. Afecta
  cuántas exposiciones caben en 30 min.

### Plan de observación (con los otros grupos y el ayudante)

- [ ] **Elegir el objeto de cada grupo** sin repetir. Los objetos que salen tarde (NGC 1097,
  1291, 1316 y 1313) van en los últimos bloques. NGC 6744 va al comienzo.
- [ ] **Elegir las dos estándares de cada noche:** una al comienzo y otra al final, para tener
  dos masas de aire distintas y poder medir el coeficiente de extinción k.
  - HD 210300 y HD 12206 obligan a exposiciones de 1–1,5 s. Mejor preferir HD 215863, HD 562 o
    HD 220881, que permiten 3–4 s.
- [ ] **Hacer la tabla del plan** con el formato del enunciado: hora, AR, Dec, objeto, tipo,
  magnitud, distancia a la Luna y observador.
- [ ] **Horario.** Confirmar con el ayudante la ventana (21:00–01:00 o hasta las 02:00) y la
  hora de inicio (el enunciado dice desde las 19:00).
- [ ] **Calibraciones.** Preguntar cuándo y cuántos bias, darks y flats se toman, y con qué
  tiempos. Tomar darks con el mismo tiempo de exposición que la ciencia o, si no, darks largos
  para escalar.
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

- [ ] **Exposición de prueba** al comienzo de cada bloque, primero en V y luego en B.
  - Revisar el **máximo** del núcleo o de la estrella y el **nivel de fondo**, ambos en ADU y
    restando el bias.
  - Si el máximo supera ~40 000 ADU, bajar el tiempo. Si el fondo es muy bajo y el máximo está
    lejos del límite, se puede subir (si el guiado lo permite).
- [ ] **Anotar en la bitácora:** tiempos de prueba, valores medidos, decisión final, hora,
  filtro, número de archivo, nubes, viento y satélites. El enunciado pide justificar los tiempos.
- [ ] **Enfoque.** Revisar que el FWHM de las estrellas no empeore durante la noche. En 2025 se
  midieron 3,9–4,4″, que es alto; puede haber influido el enfoque.
- [ ] **Satélites.** Revisar cada imagen. Con exposiciones largas, perder una cuesta más; si una
  sale contaminada, tomar otra.
- [ ] **Estándares.** Tomar 5 o más exposiciones por filtro, sin saturar y sin bajar de ~2 s si
  se puede.
- [ ] **15/10:** la Luna está alta hasta las 00:52. El fondo en B será más alto; revisarlo en la
  prueba.
- [ ] **Descargar todo** al terminar: la ciencia del grupo, las estándares y todas las
  calibraciones (bias, darks, flats y BPM).

## 3. Después, para el análisis y el informe

- [ ] **Pipeline.** Copiar `config.example.yaml` como `config.yaml` y completarlo:
  - escala de placa (0,317″/px, ya en `config.example.yaml`);
  - ganancia (0,85 e⁻/ADU);
  - objeto y coordenadas;
  - estándares con sus magnitudes B y V;
  - rutas de los archivos de las dos noches.
- [ ] **Detección con photutils.** Se usa photutils en vez de Source Extractor. El pipeline ya
  tiene las dos configuraciones (agresiva y extendida); hay que explicar en qué se diferencian
  los parámetros de deblending.
- [ ] **Astrometría.** Usar astrometry.net (modo `archivo` o `astrometry_net`); el modo `header`
  es solo aproximado.
- [ ] **Comparación de noches.** Usar `comparacion_noches.tex` (seeing, brillo del cielo,
  magnitud límite y Luna) para responder si un conjunto de imágenes es mejor que el otro.
- [ ] **Comparar lo planificado con lo obtenido:** S/N, cielo y seeing reales frente a los
  estimados en `tiempos_exposicion.md`.
- [ ] **Errores de 2025 que no hay que repetir:**
  - usar `log10` y no `log`;
  - normalizar cada imagen por su propio tiempo de exposición;
  - calcular el error del ZP como σ/√N;
  - marcar los píxeles malos como NaN, no como 0;
  - alinear las imágenes antes de combinarlas;
  - calcular el color como B − V, cruzando las fuentes por posición.
- [ ] **Código.** Va completo y comentado en el apéndice. Todos los integrantes deben poder
  explicarlo.
