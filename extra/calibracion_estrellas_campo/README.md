# Calibración con estrellas del campo (fuera de la tarea)

La noche del 15/10/2026 no hubo estrellas estándar (se observan la semana siguiente). Esta
carpeta calibra los apilados del pipeline con estrellas del propio campo de NGC 7793 que tienen
B y V de catálogo, como alternativa y para comparar después con el ZP de las estándares.

- `calibrar_campo.py`: corrige la astrometría aproximada del header contra el catálogo, hace
  fotometría de apertura en `stack_B.fits` y `stack_V.fits`, y ajusta
  `m_catálogo − m_instrumental = ZP + c·(B − V)`.
- Catálogos: APASS DR9 (VizieR II/336, por defecto) o la fotometría sintética Johnson de Gaia DR3
  (`--catalogo gaia`).

```
cd extra\calibracion_estrellas_campo
python calibrar_campo.py ..\..\pipeline\resultados\noche_15oct
python calibrar_campo.py ..\..\pipeline\resultados\noche_15oct --catalogo gaia
```

Opciones útiles: `--sin-color` (ajusta solo el ZP) y `--ra/--dec` (centro de búsqueda; por
defecto se busca el objeto por su nombre en `OBJECT`).

## Resultado: noche del 15/10/2026 (NGC 7793)

Astrometría resuelta contra el catálogo: escala 0,3163″/px y rotación −1,5°, como en 2025. Las
coordenadas `OBJCTRA`/`OBJCTDEC` del header estaban ~40′ corridas respecto de NGC 7793 (la
galaxia sí está centrada en la imagen), así que la astrometría `header` del pipeline no sirve
para esta noche.

| Filtro | Catálogo | ZP a X ≈ 1,04 (ADU/s) | Término de color c | rms | N |
|---|---|---|---|---|---|
| V | APASS | 21,469 ± 0,020 (c = 0) | −0,53 ± 0,08 si se ajusta | 0,049 | 11 |
| V | Gaia | 21,450 ± 0,018 (c = 0) | −0,13 ± 0,10 si se ajusta | 0,043 | 11 |
| B | APASS | 19,163 ± 0,096 (en B−V = 0) | **+2,01 ± 0,14** | 0,073 | 13 |
| B | Gaia | 19,239 ± 0,065 (en B−V = 0) | **+1,98 ± 0,09** | 0,078 | 11 |

- **V:** los dos catálogos coinciden a 0,02 mag. Fuera de la atmósfera (k = 0,15 supuesto) da
  21,61–21,63, 0,3 mag bajo el ZP de 2025 (21,94): pérdida de transparencia por las nubes de esa
  noche. El término de color en V no queda determinado (los catálogos discrepan); se usa c = 0.
- **B: el ZP depende fuertemente del color**, con la misma pendiente en ambos catálogos. Una
  estrella con B−V = 1,1 aparece ~2 mag más brillante en "B" de lo que corresponde a su color.
  Eso no es un término de color normal de un filtro B de Johnson (≲ 0,2): indica una **fuga al
  rojo** del filtro B (los vidrios azules dejan pasar el infrarrojo cercano, donde el CMOS es
  sensible) o que la rueda no puso el filtro B.
  - Para estrellas rojas (B−V ≈ 1,1, como las estándares K de 2025) el ZP da ~21,6 fuera de la
    atmósfera, 0,25 mag bajo 2025: lo mismo que se perdió en V por las nubes. Es decir, en 2025
    también estaba, pero no se notó porque las dos estándares eran rojas.
  - Para una estrella A0 (B−V ≈ 0) el ZP en B es ~19,4: unas 2,4 mag menos de luz que lo supuesto
    al planificar.
- **Confirmado con las estándares del 26/10** (HD 220881, B−V = 0,28, y HD 8130, B−V = 0,05):
  c_B = +2,29 ± 0,1 y ZP(B−V = 0) = 19,43 fuera de la atmósfera, de acuerdo con este ajuste
  (19,42–19,50). En V, el ZP de las estándares es 21,75: el 15/10 perdió ~0,13 mag por nubes, no
  0,3 como sugería la comparación con 2025. El pipeline de la tarea ya aplica el término de color
  (`termino_color` en `config.yaml`). Falta preguntar en el observatorio qué filtro B tiene la rueda.

### Catálogo calibrado y diagrama color-magnitud

`calibrar_catalogos.py` calibra el catálogo "agresiva" del pipeline (B y V cruzados por
posición, S/N > 5) y dibuja el diagrama color-magnitud:

```
python calibrar_catalogos.py ..\..\pipeline\resultados\noche_15oct
```

- El ZP se ajusta sobre las magnitudes instrumentales del propio catálogo, con las mismas 13
  estrellas de APASS (el catálogo usa otra apertura y otra corrección de apertura que
  `calibrar_campo.py`): V = 21,544 ± 0,028 (sin término de color, rms 0,075) y B = 19,316 ±
  0,114 con c_B = +2,03 ± 0,15 (rms 0,10).
- Color: `B−V = [(b − v) + ZP_B − ZP_V] / (1 − c_B)`. Con c_B ≈ 2 el denominador es ≈ −1: el
  color instrumental b − v **disminuye** cuando la estrella es más roja, por la fuga al rojo.
- Comprobación con las 13 estrellas de APASS: ΔV = −0,01 (rms 0,14) y Δ(B−V) = −0,02 (rms 0,10).
- 97 fuentes con B y V: 83 del campo y 14 dentro del disco de NGC 7793 (elipse de 10′ × 6′, AP
  84°). A la distancia de NGC 7793 (~3,9 Mpc) las fuentes del disco son cúmulos, regiones HII o
  mezclas de estrellas, no estrellas individuales.
- **Validez:** el color solo está calibrado para B−V ≈ 0,45–1,12 (rango de las estrellas de
  APASS; franja gris del diagrama). Los colores fuera de ese rango, y en particular los muy rojos
  (B−V > 1,5) de las fuentes débiles, son extrapolaciones de un término de color anómalo y no
  deben interpretarse como colores físicos hasta confirmar el filtro B.

Salidas: `resultados/catalogo_calibrado_noche_15oct.csv` y
`resultados/cmd_ngc7793_noche_15oct.png`.

### Astrometría: astrometry.net frente al ajuste con Gaia

`comparar_astrometria.py` compara la solución de astrometry.net (que el pipeline guardó en el
header del apilado) con la que se ajusta aquí contra Gaia DR3:

| | Escala (″/px) | Rotación | Centro de la imagen | Residuo contra Gaia (41 estrellas) |
|---|---|---|---|---|
| astrometry.net | 0,31616 | −1,50° | 23h57m45,43s −32°34′59,7″ | 0,75″ rms |
| Ajuste con Gaia | 0,31631 | −1,52° | 23h57m45,43s −32°34′59,8″ | 0,31″ rms |

Las dos soluciones coinciden: 0,3″ de diferencia mediana en la imagen (3,4″ como máximo en las
esquinas, donde astrometry.net agrega distorsión SIP). La escala coincide con la medida en 2023 y
2025 (0,3168 y 0,3163″/px). El centro de NGC 7793 queda en el píxel (2223, 1964).

Las tablas por estrella (`estrellas_campo_*.csv`) traen `dm_B` y `dm_V` (catálogo −
instrumental) para cada estrella usada. El ZP sale a la masa de aire del apilado; para pasarlo a
fuera de la atmósfera se usa un k supuesto (B 0,25, V 0,15), porque sin estándares a dos masas
de aire no se puede medir.
