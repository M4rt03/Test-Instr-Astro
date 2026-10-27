# Observaciones: noche del 15/10/2026 (resumen para el informe)

**Condiciones y datos.** La noche del 8/10 se perdió por nubes. El 15/10 se observó con nubes
intermitentes y viento, con la Luna al 26 % a 81° del objeto. De NGC 7793 se obtuvieron 5
exposiciones de 120 s en B y 4 de 120 s más 2 de 60 s en V (22:41–23:06, hora local, masa de
aire 1,03–1,05). NGC 1291 quedó con una sola exposición casi sin señal: se nubló y el
observatorio cerró el domo. No se observaron estrellas estándar ni se tomaron bias; como nivel
cero se usó el master dark de 120 s (que incluye el bias) y los flats de crepúsculo los tomó el
ayudante. Se excluyeron tres exposiciones: una por pérdida de transparencia (0,28 mag) y dos con
estrellas dobles por movimientos del telescopio con el viento. Los apilados finales suman 480 s
en B (4 × 120 s) y 360 s en V (3 × 120 s). Las estrellas estándar se observaron el 26/10 (ver
abajo).

**Seeing.** El FWHM medido en los apilados es 3,7″ en B y 4,1″ en V (3,2–4,1″ por exposición),
con estrellas alargadas (elipticidad mediana 0,3–0,4). El observatorio reportó un seeing de
1,9–2,1″ para el sitio, así que la diferencia se atribuye al viento sobre el telescopio y no a la
atmósfera. El límite de difracción del telescopio (0,19″ en B y 0,23″ en V) es despreciable
frente a estos valores.

**Cielo.** El fondo fue de 0,14 ADU/s/píxel en ambos filtros, que corresponde a
μV ≈ 21,1–21,3 mag/arcsec² en V (según se use el zeropoint de las estrellas del campo o el de las
estándares). En B no se da un valor en magnitudes porque el zeropoint en B depende fuertemente
del color (ver abajo).

**Astrometría.** Las coordenadas del header (`OBJCTRA`/`OBJCTDEC`) quedaron ~40′ desplazadas de
NGC 7793, aunque la galaxia está centrada en la imagen, así que se resolvió la astrometría con
astrometry.net: escala 0,3162″/px, rotación −1,50°, centro 23h57m45,4s −32°35′00″. Un ajuste
independiente contra 41 estrellas de Gaia DR3 da 0,3163″/px y −1,52°, con 0,3″ de diferencia
mediana entre ambas soluciones; la escala coincide con la medida en 2023 y 2025.

**Estrellas estándar (26/10).** Se observaron HD 220881 (A9V, B−V = 0,28; 5 × 40 s en B y
5 × 8 s en V) y HD 8130 (A0V, B−V = 0,05; 5 × 70 s en B y 7 × 5 s en V), con darks del mismo
tiempo de exposición y los flats del 15/10, sin saturar (picos de 16 000–33 000 ADU). Con las
magnitudes Johnson de Mermilliod y k_B = 0,25, k_V = 0,15 (valores típicos: cada estrella se
observó una sola vez y la extinción no se pudo medir):

- **V:** ZP = 21,750 ± 0,012 ADU/s fuera de la atmósfera; las dos estrellas coinciden a 0,02 mag.
- **B:** la estrella más azul da un ZP 0,54 mag menor. Ajustando m_B − b = ZP + c (B−V), resulta
  ZP(B−V = 0) = 19,43 ± 0,01 y c = +2,29 ± 0,06 (± 0,1 si se incluye la incertidumbre de k_B).

El pipeline aplica este zeropoint y el término de color a los catálogos del 15/10. Comparado con
13 estrellas de APASS del campo, el B−V calibrado coincide a −0,04 mag (rms 0,07), aunque el
término de color se midió con estrellas de B−V ≤ 0,28; las magnitudes V y B salen ~0,15 mag más
débiles que APASS, que es la pérdida de transparencia por las nubes del 15/10 (el zeropoint viene
de una noche despejada).

**Zeropoint con estrellas del campo.** Antes de tener estándares, el zeropoint se estimó con 10–13
estrellas del campo (11 < V < 16, fuera de la galaxia) con magnitudes de APASS DR9 y de la
fotometría sintética Johnson de Gaia DR3. En V ambos catálogos coinciden: ZP = 21,46 ± 0,02
ADU/s a masa de aire 1,04 (rms 0,04–0,05 mag), que fuera de la atmósfera (k_V = 0,15 supuesto)
es 21,62: 0,13 mag bajo el de las estándares del 26/10 (21,75), por las nubes.

**Filtro B.** En B la diferencia entre la magnitud de catálogo y la instrumental crece con el
color de la estrella, con un término de color de +2,0 ± 0,1 tanto con APASS como con Gaia (rms
0,07 mag). Un filtro B de Johnson tiene términos de color de ≲ 0,2, por lo que el resultado
indica una fuga hacia el rojo del filtro B (o un filtro distinto en la rueda): las estrellas
rojas aparecen mucho más brillantes en "B" de lo que corresponde. Para una estrella roja
(B−V ≈ 1,1, como las estándares de 2025) el zeropoint resultante coincide con el de 2025
descontando las nubes, lo que sugiere que el efecto ya existía y no se notó porque esas
estándares eran rojas. Las estándares del 26/10 lo confirman (c = +2,29): con ese término de
color, las dos estándares de 2025 (B−V ≈ 1,07 y 0,96) debían dar ZP de 21,89 y 21,63, y dieron
21,99 y 21,70. La diferencia de 0,29 mag entre ellas, que en 2025 se atribuyó a una serie mal
tomada, es el término de color. Falta preguntar en el observatorio qué filtro hay en la posición
B de la rueda.
