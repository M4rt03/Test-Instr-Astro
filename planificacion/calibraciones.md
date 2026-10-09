# Calibraciones y estrellas estándar: tiempos para las noches del 15 y 22/10

Plan de bias, darks, flats y estrellas estándar para el MAS500 (El Sauce), con los objetos del
grupo: NGC 7793 y NGC 1291. Todas las horas son de Chile (UTC−3).

La noche del 8/10 se nubló: quedan el 15/10 y, como respaldo, el 22/10. Los flats los toma el
ayudante.

## Resumen

| Imagen | Cuántas | t_exp | Cuándo | Tiempo aprox. |
|---|---|---|---|---|
| Bias | 25 al inicio y 25 al final | 0 s | Al comienzo y al terminar | < 1 min cada serie (en 2025, un bias cada 0,7 s) |
| Darks | 10 por cada t_exp de ciencia: 120, 90 y 60 s | 120 / 90 / 60 s | Con la cúpula cerrada: antes de anochecer o al final | ~45 min (mínimo: 10 × 120 s, ~21 min) |
| Flats de cielo | 10–15 por filtro, primero B y luego V | 1–30 s, para 20 000–30 000 ADU | **Los toma el ayudante.** En el crepúsculo serían 20:05–20:40 (15/10) y 20:10–20:45 (22/10) | ~35 min |
| Estándares | 5 por filtro en cada visita, 5 visitas | 1–3 s (tabla abajo) | Comienzo, mitad y final de la noche (tabla de visitas) | ~25 min en total |

Las calibraciones sirven a todos los grupos de la misma noche. Conviene repartirlas: bias y flats
una sola vez, y los darks con los t_exp de todos los grupos.

## Horario del Sol en El Sauce

| Evento | 15/10 | 22/10 |
|---|---|---|
| Puesta del Sol | 19:54 | 19:59 |
| Sol a −6° (fin del crepúsculo civil) | 20:19 | 20:24 |
| Sol a −12° (fin del náutico) | 20:47 | 20:53 |
| Sol a −18° (comienza la noche astronómica) | 21:17 | 21:24 |

Calculado para −30,4597°, −70,7503° (header de 2025). Coincide con Staralt para Tololo: el 8/10
la noche empezaba a las 0h09m UT, es decir, a las 21:09.

### Secuencia de la noche del 15/10

| Hora | Qué |
|---|---|
| Antes de 20:00 | Enfriar la cámara y esperar que `CCD-TEMP` se estabilice en la consigna. Bias (25). Darks, si se pueden tomar con la cúpula cerrada |
| 20:05–20:45 | Flats (los toma el ayudante) |
| 20:45–21:15 | Enfoque y apuntado; exposiciones de prueba |
| 21:15–21:35 | **Estándares, visita 1:** HD 202941 y HD 8130 |
| 22:30–23:30 | NGC 7793 |
| ~23:30 | **Estándares, visita 2:** HD 220881 |
| 23:50–00:45 | NGC 1291 (o 01:00–02:00 si la ventana llega a las 02:00) |
| ~00:45–01:00 | **Estándares, visita 3:** HD 8130 y HD 202941 |
| Al terminar | Bias (25) y los darks que falten |

La Luna (26 %) se pone a la 01:00 y está a más de 79° de los dos objetos: casi no afecta.

### Secuencia de la noche del 22/10 (respaldo)

Todo ocurre ~30 min antes que el 15/10. La Luna (88 %) está sobre 38° toda la ventana; ver la
S/N esperada en [`tiempos_exposicion.md`](tiempos_exposicion.md) (sección 4).

| Hora | Qué |
|---|---|
| Antes de 20:05 | Cámara a la consigna y bias (25) |
| 20:10–20:50 | Flats (los toma el ayudante) |
| 20:50–21:25 | Enfoque y apuntado; exposiciones de prueba |
| 21:25–21:45 | **Estándares, visita 1:** HD 202941 y HD 8130 |
| 21:50–22:50 | NGC 7793 (está a ~31° de la Luna toda la noche; la hora casi no cambia la S/N) |
| ~23:00 | **Estándares, visita 2:** HD 220881 |
| 23:45–00:40 | NGC 1291 (mejor 01:00–02:00 si la ventana llega a las 02:00) |
| ~00:45 | **Estándares, visita 3:** HD 8130 y HD 202941 |
| Al terminar | Bias (25) y los darks que falten |

## 1. Bias

- **25 exposiciones de 0 s** al comienzo y 25 al final. Con 25, el ruido del master bias es
  3,56 e⁻ / √25 = 0,7 e⁻, despreciable frente al ruido de lectura de cada imagen.
- Comparar el nivel del bias del inicio con el del final. En 2025 fue 92,5 ADU; si cambia más de
  ~1 ADU durante la noche, usar el bias más cercano a cada imagen.

## 2. Darks

- **Mismo t_exp que la ciencia:** 120 s (NGC 7793 en B), 90 s (NGC 7793 en V y NGC 1291 en B) y
  60 s (NGC 1291 en V). Tomar 10 de cada uno.
- **Misma temperatura del sensor que la ciencia.** En 2025 la cámara quedó a −14 °C sin llegar a la
  consigna de −16 °C. Revisar `CCD-TEMP` en el header.
- **Por qué con el mismo t_exp:** en 2025 los darks de 100 s quedaron 3,9 ADU *bajo* el bias. En
  este sensor CMOS, bias + corriente oscura escalada no reproduce bien un dark largo. Con darks
  del mismo t_exp, lo que se resta es exactamente lo que tiene la imagen.
- **Cuándo:** son ~45 min, así que no conviene tomarlos de noche. Se pueden tomar con la cúpula
  cerrada y el telescopio tapado, antes de anochecer o al terminar, siempre a la misma
  temperatura.
- **Mínimo, si no hay tiempo:** 10 darks de 120 s. El pipeline (`make_dark_current`) calcula
  (dark − bias)/t y lo escala a cada t_exp. Así el error en las imágenes de 60–90 s queda bajo
  ~1 ADU, y el fondo local se resta igual en la fotometría.
- **Estándares y flats (1–30 s):** no necesitan darks propios; la corriente oscura en ese tiempo
  es < 1 e⁻.

## 3. Flats de cielo (crepúsculo)

**Los toma el ayudante.** Lo que sigue sirve para pedirle lo necesario y revisar que los flats
cumplan: B y V, 10–15 por filtro, con 20 000–30 000 ADU, y saber de qué noche son.

**Nivel:** 20 000–30 000 ADU sobre el bias, dentro del rango lineal (límite de trabajo: 40 000 ADU).
- Con 25 000 ADU, cada flat tiene 21 000 e⁻ por píxel: 0,7 % de ruido.
- Al combinar 10–15 flats, el ruido baja a ~0,2 %, menos que la precisión de la fotometría.
- Descartar los que queden bajo 10 000 ADU o sobre 40 000 ADU. El pipeline avisa cuando un flat
  supera el 80 % de la saturación.

**Cuándo:**
- Con 0,317″/px, cada píxel cubre poco cielo. Para 25 000 ADU en 5 s hace falta un cielo de
  ~10 mag/arcsec² (con el ZP medido).
- Eso ocurre aproximadamente con el Sol entre −3° y −8°: **20:00–20:35 el 8/10** y
  **20:05–20:40 el 15/10**. Es una estimación: la exposición de prueba manda.
- **Primero B y luego V.** B recibe menos luz (eficiencia 0,22 contra 0,37 en V), así que se
  aprovecha cuando el cielo está más brillante.

**Procedimiento en TheSkyX:**
1. Apuntar ~15° del cenit hacia el este (lado opuesto al Sol). Es la zona más uniforme del
   cielo de crepúsculo.
2. Tomar una prueba de 1 s en B en cuanto el cielo deje de saturar.
3. Calcular el siguiente t_exp con `t_nuevo = t × 25 000 / nivel`. Agregar ~20 % por cada minuto
   que pase, porque el cielo se oscurece rápido: ~0,2 mag por minuto con el Sol bajando ~0,2° por
   minuto.
4. Mover el telescopio ~1′ entre exposiciones (*dithering*), para que la mediana elimine las
   estrellas.
5. Con 10–15 flats en B (o si B pasa de ~30 s), cambiar a V y repetir, partiendo de nuevo con
   una prueba.
6. Anotar el nivel de cada flat en la bitácora.

**Alternativa:** si el cielo está nublado al atardecer, preguntar si hay pantalla o lámpara para
flats de cúpula. Si no, se usan los flats de la otra noche: el filtro y la óptica no cambian.

## 4. Estrellas estándar

### Catálogo: astroquery + Mermilliod (Johnson UBV)

[`estrellas_estandar.py`](estrellas_estandar.py) consulta VizieR con astroquery y elige, para cada
estrella, B y V en el sistema de Johnson:

1. **II/168, Mermilliod (1991):** medias homogéneas de fotometría fotoeléctrica UBV de Johnson.
   Es el catálogo estándar para estrellas brillantes.
2. **I/239, Hipparcos,** si la estrella no está en Mermilliod. Hipparcos indica el origen de
   cada valor: G = fotometría terrestre (Johnson), H = V calculada desde la banda Hp, T = B−V
   calculado desde Tycho.

Catálogos descartados:
- **Tycho (BT, VT):** no es Johnson. Usarlo fue uno de los errores de 2025. El pipeline puede
  convertirlo, pero con ~0,02–0,05 mag de error.
- **APASS (II/336):** tiene B y V, pero satura con estas estrellas. Difiere hasta 0,43 mag:

  | Estrella | V Mermilliod / Hipparcos | V APASS | B APASS |
  |---|---|---|---|
  | HD 210300 | 6,440 | 6,874 | 6,929 (Merm.: 6,590) |
  | HD 202941 | 7,070 | 7,170 | 7,069 |
  | HD 12206 | 6,790 | 6,882 | 6,802 |
  | HD 225200 | 6,380 | 7,932 | 7,841 |

### Resultado

| Estrella | Uso | Tipo | Fuente | Calidad | B | V | t_exp B | t_exp V |
|---|---|---|---|---|---|---|---|---|
| HD 202941 | seleccionada | A0V | Hipparcos (terrestre) | alta | 7.072 | 7.070 | 1.5 s | 1.0 s |
| HD 207480 | seleccionada | A1V | Mermilliod (7 obs.) | alta | 7.190 | 7.140 | 1.5 s | 1.5 s |
| HD 210300 | seleccionada | A5V | Mermilliod (4 obs.) | alta | 6.590 | 6.440 | 1.0 s | 1.0 s |
| HD 212643 | seleccionada | A0V | Mermilliod (4 obs.) | alta | 6.260 | 6.290 | 1.0 s | 1.0 s |
| HD 220881 | seleccionada | A9V | Mermilliod (4 obs.) | alta | 7.720 | 7.440 | 3.0 s | 1.5 s |
| HD 562 | seleccionada | A2V | Mermilliod (8 obs.) | alta | 7.787 | 7.650 | 3.0 s | 2.0 s |
| HD 8130 | seleccionada | A0V | Mermilliod (9 obs.) | alta | 7.493 | 7.448 | 2.5 s | 2.0 s |
| HD 12206 | seleccionada | A0V | Mermilliod (4 obs.) | alta | 6.810 | 6.790 | 1.0 s | 1.0 s |
| HD 223884 | reserva | A5V | Mermilliod (4 obs.) | alta | 6.430 | 6.240 | 1.0 s | 1.0 s |
| HD 225200 | reserva | A1V | Hipparcos (V de Hp, B−V terrestre) | media | 6.386 | 6.380 | 1.0 s | 1.0 s |
| HD 7323 | reserva | A0V | Mermilliod (7 obs.) | alta | 7.928 | 7.830 | 3.5 s | 2.5 s |
| HD 195500 | solo extinción | A1V | Hipparcos (V de Hp, B−V de Tycho) | **baja** | 7.380 | 7.320 | 2.0 s | 1.5 s |
| HD 215863 | solo extinción | A2V | Hipparcos (V de Hp, B−V de Tycho) | **baja** | 7.836 | 7.690 | 3.0 s | 2.5 s |

- En las estrellas que están en ambos catálogos, Mermilliod e Hipparcos coinciden a 0,01 mag.
- **HD 195500 y HD 215863** no tienen fotometría Johnson medida: su B−V viene de Tycho. Por eso
  se **reemplazaron** por estrellas de Mermilliod en la misma zona del cielo:
  - **HD 207480** (AR 21:49) en lugar de HD 195500;
  - **HD 212643** (AR 22:26) en lugar de HD 215863. Es brillante: en 1 s llega a ~30 000–35 000
    ADU con seeing de 2″, así que hay que revisar el máximo en la prueba.

  Las reemplazadas igual sirven para medir la extinción, porque ahí se comparan magnitudes de la *misma*
  estrella y su magnitud de catálogo se cancela.
- **t_exp:** es ~60 % del tiempo de saturación con seeing de 2,0″ (mismo modelo de
  [`tiempos_exposicion.md`](tiempos_exposicion.md)). Tomar 5 o más por filtro, y más si son de 1 s.

### Visitas: dos masas de aire para medir la extinción

El pipeline usa un k fijo (`extincion` en `config.yaml`: 0,25 en B y 0,15 en V). Si se observa la
misma estrella a dos masas de aire, k se puede medir:

`k = Δm_instrumental / ΔX`

Basta con dos estrellas, cada una al comienzo y al final:

| Visita | Hora 15/10 | Hora 22/10 | Estrellas | X el 15/10 | X el 22/10 |
|---|---|---|---|---|---|
| 1 | 21:15–21:35 | 21:25–21:45 | HD 202941 · HD 8130 | 1,00 · 1,61 | 1,01 · 1,37 |
| 2 | ~23:30 | ~23:00 | HD 220881 (o HD 562) | 1,00 | 1,00 |
| 3 | ~00:45–01:00 | ~00:45 | HD 8130 · HD 202941 | 1,01 · 1,46 | 1,01 · 1,54 |

- HD 8130 pasa de X ≈ 1,4–1,6 a 1,0 y HD 202941 de 1,0 a ~1,5: ΔX entre 0,4 y 0,6.
- Si la ventana llega a las 02:00, dejar la visita 3 para el final: el ΔX sube (HD 202941 a 2,0–2,4).
- La Luna no afecta a las estándares: son de V 6–8 y las exposiciones son de 1–3 s.
- Cada visita toma ~4–5 min por estrella: apuntar, centrar y tomar 5 + 5 exposiciones.

## 5. Cómo correr el script

En Anaconda Prompt (una sola vez):

```
pip install astroquery
```

Luego:

```
cd C:\Users\marti\Test-Instr-Astro\planificacion
python estrellas_estandar.py
python estrellas_estandar.py --todas --csv estandares_bv.csv
```

Imprime la tabla de arriba y el bloque `estandares:` listo para `pipeline/config.yaml`, con
coordenadas, magnitudes B y V, y rutas de archivos. Cambiar `--noche noche2` para la segunda noche.

## Falta confirmar con el ayudante

- Si los darks se pueden tomar con la cúpula cerrada antes de anochecer, y si hay una biblioteca
  de darks a la temperatura de consigna.
- Qué noche toma los flats y con qué nivel (20 000–30 000 ADU), para usarlos con las imágenes
  correctas.
- Qué t_exp usan los otros grupos, para tomar sus darks en la misma serie.
