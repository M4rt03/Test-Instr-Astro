# Observaciones: noche del 15/10/2026 (resumen para el informe)

**Condiciones y datos.** La noche del 8/10 se perdió por nubes. El 15/10 se observó con nubes
intermitentes y viento, con la Luna al 26 % a 81° del objeto. De NGC 7793 se obtuvieron 5
exposiciones de 120 s en B y 4 de 120 s más 2 de 60 s en V (22:41–23:06, hora local, masa de
aire 1,03–1,05). NGC 1291 quedó con una sola exposición casi sin señal: se nubló y el
observatorio cerró el domo. No se observaron estrellas estándar ni se tomaron bias; como nivel
cero se usó el master dark de 120 s (que incluye el bias) y los flats de crepúsculo los tomó el
ayudante. Se excluyeron tres exposiciones: una por pérdida de transparencia (0,28 mag) y dos con
estrellas dobles por movimientos del telescopio con el viento. Los apilados finales suman 480 s
en B y 420 s en V.

**Seeing.** El FWHM medido en los apilados es 3,7″ en B y 4,1″ en V (3,2–4,1″ por exposición),
con estrellas alargadas (elipticidad mediana 0,3–0,4). El observatorio reportó un seeing de
1,9–2,1″ para el sitio, así que la diferencia se atribuye al viento sobre el telescopio y no a la
atmósfera. El límite de difracción del telescopio (0,19″ en B y 0,23″ en V) es despreciable
frente a estos valores.

**Cielo.** El fondo fue de 0,14 ADU/s/píxel en ambos filtros. Con el zeropoint de la noche
corresponde a μV ≈ 21,1 mag/arcsec² en V. En B no se da un valor en magnitudes porque el
zeropoint en B depende fuertemente del color (ver abajo).

**Astrometría.** Las coordenadas del header (`OBJCTRA`/`OBJCTDEC`) quedaron ~40′ desplazadas de
NGC 7793, aunque la galaxia está centrada en la imagen, así que se resolvió la astrometría con
astrometry.net: escala 0,3162″/px, rotación −1,50°, centro 23h57m45,4s −32°35′00″. Un ajuste
independiente contra 41 estrellas de Gaia DR3 da 0,3163″/px y −1,52°, con 0,3″ de diferencia
mediana entre ambas soluciones; la escala coincide con la medida en 2023 y 2025.

**Zeropoint con estrellas del campo.** Sin estándares, el zeropoint se estimó con 10–13
estrellas del campo (11 < V < 16, fuera de la galaxia) con magnitudes de APASS DR9 y de la
fotometría sintética Johnson de Gaia DR3. En V ambos catálogos coinciden: ZP = 21,46 ± 0,02
ADU/s a masa de aire 1,04 (rms 0,04–0,05 mag), que fuera de la atmósfera (k_V = 0,15 supuesto)
es 21,62, unas 0,3 mag bajo el valor medido con estándares en 2025 (21,94); la diferencia es
consistente con la pérdida de transparencia por las nubes.

**Filtro B.** En B la diferencia entre la magnitud de catálogo y la instrumental crece con el
color de la estrella, con un término de color de +2,0 ± 0,1 tanto con APASS como con Gaia (rms
0,07 mag). Un filtro B de Johnson tiene términos de color de ≲ 0,2, por lo que el resultado
indica una fuga hacia el rojo del filtro B (o un filtro distinto en la rueda): las estrellas
rojas aparecen mucho más brillantes en "B" de lo que corresponde. Para una estrella roja
(B−V ≈ 1,1, como las estándares de 2025) el zeropoint resultante coincide con el de 2025
descontando las nubes, lo que sugiere que el efecto ya existía y no se notó porque esas
estándares eran rojas. Mientras no se confirme con estándares de colores distintos, las
magnitudes B y los colores B−V de esta noche deben considerarse no calibrados.
