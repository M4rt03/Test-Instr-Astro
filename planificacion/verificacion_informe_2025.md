# Verificación de los datos del instrumento y del informe de 2025

Revisión de los datos de la cámara y de los cálculos del informe "Procesamiento y Calibración
Astrométrica y Fotométrica de NGC 6618" (2025), hecha el 7/10/2026.

**Fuentes usadas:**
- Los headers FITS que quedaron impresos en los notebooks del repositorio: datos de 2025
  (NGC 6618) y de 2023 (HIP 117445, NGC 2070, NGC 6514).
- Las tablas de fuentes de esos notebooks.
- Gaia DR3 y SIMBAD (vía CDS).
- Cálculo propio de posiciones del Sol y la Luna (fórmulas de Meeus).

## 1. Cámara y telescopio

### Confirmado en los headers (escritos por TheSkyX)

| Dato | Valor | Keyword |
|---|---|---|
| Tamaño de imagen | 4096 × 4096 px | `NAXIS1`, `NAXIS2` |
| Tamaño de píxel | 9 µm | `XPIXSZ`, `YPIXSZ` |
| Binning | 1 × 1 | `XBINNING` |
| Ganancia | 0,85 e⁻/ADU | `GAIN` |
| Profundidad | 16 bits (0–65 535 ADU) | `BITPIX = 16`, `BZERO = 32768` |
| Temperatura del sensor | −14,3 °C con consigna −16 °C (2025); −24,8 °C con consigna −25 °C (2023) | `CCD-TEMP`, `SET-TEMP` |
| Apertura | 500 mm | `APTDIA` |
| Focal | 6500 mm (valor ingresado en el software, no medido) | `FOCALLEN` |
| Sitio (El Sauce) | −30,4597°, −70,7503°, 1600 m | `OBSGEO-B/L/H` |
| Hora | `DATE-OBS` en **UTC** ("UTC of start exp.") | `DATE-OBS` |
| Software | TheSkyX 10.5.0 | `SWCREATE` |

### Confirmado en la ficha de Moravian

Fuente: [C4 Series CMOS Cameras](https://www.gxccd.com/art?id=607&lang=409), consultada el 7/10/2026.

| Dato | Valor |
|---|---|
| Sensor | Gpixel GSENSE4040, 4096 × 4096 px de 9 µm (~37 × 37 mm). Coincide con el header |
| Versión | Dos versiones: FSI y BSI. La ganancia de 0,85 e⁻/ADU del header corresponde a la **FSI** (la BSI tiene 0,37) |
| Modo de lectura | **16-bit HDR**: combina los canales de 12 bits de baja y alta ganancia; la ganancia no se puede elegir |
| Capacidad de pozo (FSI, HDR) | 56 600 e⁻, unos 66 600 ADU: el ADC (65 535) satura antes que el pozo |
| Ruido de lectura (FSI, HDR) | 3,9 e⁻ RMS |
| Otros modos (FSI) | 12 bits alta ganancia: 0,85 e⁻/ADU, 3,9 e⁻, 3540 e⁻ · 12 bits baja ganancia: 19,5 e⁻/ADU, 34,5 e⁻, 80 000 e⁻ |
| Obturador | Electrónico tipo *rolling shutter*: cada fila empieza 21 µs después que la anterior (84,5 ms entre la primera y la última); todas se exponen el mismo tiempo. El obturador mecánico solo se usa para darks y bias |
| Exposición mínima | 21 µs |
| Descarga | 0,25 s (USB 3) o 1,6 s (USB 2) en modo HDR |
| Enfriamiento | Regulado hasta 33 °C bajo el ambiente (versión EC) o 28 °C (estándar), con precisión de 0,1 °C |
| Binning | Por software (el de hardware, 2 × 2, da peor calidad) |
| Píxeles defectuosos | Sensor grado 1: hasta 300 píxeles defectuosos |

### Sin confirmar

La ficha no publica la **eficiencia cuántica** ni la **corriente oscura** del C4-16000. Siguen
como supuestos en `tiempos_exposicion.py`:
- la corriente oscura se puede medir con los darks de 2025 (`paso_calibracion` → `calibracion.json`);
- la eficiencia cuántica queda incluida en el zeropoint, que se mide con las estándares.

## 2. Escala de placa: tres valores distintos

| Origen | Escala | Campo (4096 px) |
|---|---|---|
| Header: 206 265 × 9 µm / 6500 mm | 0,2856″/px | 19,5′ |
| Informe 2025 (`PIXSCALE = 0.36`, escrito a mano en su código) | 0,36″/px | 24,6′ |
| Matriz CD de astrometry.net en el informe (Tabla 4) | 0,3163″/px | 21,6′ |
| **Medición propia con Gaia ([`medir_escala.py`](medir_escala.py))** | **0,3168″/px** | **21,6′** |

### Cómo se midió

1. **Datos:** `Calibracao.ipynb` tiene impresa la tabla completa de las 8 fuentes detectadas en
   una exposición de 9 s en B de la estándar HIP 117445 (2023-10-06). La fuente más grande, en
   (2042, 2036), es la estrella, en el centro de la imagen.
2. **Catálogo:** se descargaron de VizieR las estrellas de Gaia DR3 con G < 15 del campo y se
   corrigió su movimiento propio a la fecha de la imagen.
3. **Búsqueda:** se probó con escalas de 0,24 a 0,40″/px, rotaciones de 0 a 360° y ambas
   orientaciones, contando cuántas de las otras 7 fuentes caen a menos de 6″ de una estrella
   de Gaia. Solo una solución hace coincidir las 7:

   | Escala | Estrellas que coinciden (de 7) |
   |---|---|
   | 0,2856″/px | 1 |
   | **0,3168″/px** | **7** |
   | 0,36″/px | 1–2 |

4. **Ajuste final:** escala 0,3168″/px (x: 0,3165, y: 0,3171), rotación −0,5° y residuo rms de
   0,84″, consistente con centroides de segmentos de 5–15 píxeles.

### Conclusión

La escala es **0,317″/px** y el campo, **21,6′ × 21,6′**.
- La medición de 2023 y la solución de astrometry.net de 2025 coinciden, así que la óptica no
  cambió entre esos años.
- El `FOCALLEN = 6500` del header es nominal: la focal efectiva es ~5860 mm (f/11,7),
  probablemente por un corrector o reductor.
- El 0,36″/px del informe de 2025 no tiene respaldo y sobreestima el campo en un 14 %.

## 3. Revisión de los cálculos del informe

| Punto del informe | Valor del informe | Valor verificado | Estado |
|---|---|---|---|
| Campo de visión | 4096 × 0,36″ = 24,58′ | La aritmética es correcta, pero la escala no: 21,6′ | ✗ |
| Escala (Tabla 4, matriz CD) | No la calculan | 0,3163″/px, rotación 1,3° | Dato correcto, no aprovechado |
| Seeing HIP 117678 B | 13,17 px → 4,47″ | 13,17 × 0,36 = 4,74″ (error de transcripción) | ✗ |
| Seeing HIP 117678 V, HIP 116375 B y V | 12,42 → 4,47″; 10,73 → 3,86″; 12,42 → 4,47″ | Aritmética correcta con 0,36 | ✓ |
| Seeing promedio | 4,39″ | 4,39″ con 0,36 (usa el 4,74 correcto). Con 0,317″/px: **3,86″** | Escala errada |
| Horas de las tablas 1 y 2 | Rotuladas como CLT (UTC−3) | Son **UTC**. M17 se observó a las 20:46–20:57 hora local; HIP 116375 a las 22:08 y HIP 117678 a las 23:52–23:54 | ✗ |
| Culminación de M17 (2/10) | ~22:30, 70° | 19:18 hora local, 75,7° | ✗ |
| Altura de M17 durante la observación | No la dan | 64–65° (masa de aire 1,10–1,11), igual al `CENTALT` y `AIRMASS` del header. Si la hora fuera 23:46 local, M17 habría estado a 27° | — |
| Culminación 24/9 | Ambos objetos ~22:30, 70–80° | NGC 7009: 22:32, 71° ✓. NGC 7293: 23:57, 80° | Parcial |
| Puesta de Sol | 19:44 (24/9), 19:50 (2/10) | 19:41 y 19:46 en El Sauce | ≈ ✓ |
| Luna el 2/10 | 78,8 % | 78,8 %, a 41° de M17 | ✓ |
| HIP 116375: V, M, d | 7,08; 0,84; 177,7 pc = 579,5 ly | V = 7,077 (de Tycho); M = 0,83; 579,6 ly | ✓ |
| HIP 117678: V, M, d | 7,06; 1,17; 150,9 pc = 492,0 ly | V = 7,058; M = 1,17; 492,2 ly | ✓ |
| B−V (Johnson) | 1,075 y 0,972 | 0,85·(BT−VT) = 1,073 y 0,961 | ✓ / ≈ |
| Magnitudes usadas para el ZP | Tycho (BT, VT) | Debían ser Johnson; además, las etiquetas de las dos estrellas están cruzadas en el código | ✗ |
| Zeropoints | B 39,26 · V 39,66 | Imposibles: una estrella de V ≈ 7 daría 10¹³ ADU/s. Ni con los 16,7 millones de píxeles saturados en 3 s se llega a 4·10¹¹ ADU/s | ✗ |
| Magnitudes del catálogo | `np.log` | Debía ser `np.log10` | ✗ |
| Error del ZP | σ/N | σ/√N | ✗ |
| Bias | "exposición de 1 s" en el texto, 0 s en la Tabla 3 | Inconsistente | Menor |
| Tipo de detector | "cámara CCD" en varias partes | Es CMOS | Menor |
| Darks | Tomados el 25/9, una semana antes de la ciencia | En 2025 el sensor no llegó a la consigna (−14,3 °C frente a −16 °C), así que la corriente oscura puede diferir | Revisar |

**Límite de difracción** (no lo calcula el informe): 1,22 λ/D = 0,22″ en B (440 nm) y 0,28″ en V
(550 nm). El seeing medido (~3,9″) es unas 14–18 veces mayor, con un muestreo de ~12 píxeles por FWHM.

## 4. Consecuencias para 2026

- **Escala y campo:** usar 0,317″/px; ya está en `pipeline/config.example.yaml` y en
  `tiempos_exposicion.py`. El campo de 21,6′ hace muy justas a NGC 247, NGC 300 y la Hélice.
- **Hora:** la de los FITS está en UTC; la hora local se obtiene restando 3 h. El pipeline ya
  lo trata así.
- **Sitio:** las coordenadas del header (−30,4597°, −70,7503°) quedan en `config.example.yaml`.
  Para la visibilidad, el enunciado pide Cerro Tololo.
- **Ruido de lectura, corriente oscura y zeropoint:** siguen pendientes; medirlos con los datos
  de 2025 (ver [`pendientes.md`](pendientes.md)).
