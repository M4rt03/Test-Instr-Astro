# FIS2421 – Tarea 2: candidatos a observar

**Sitio:** Cerro Tololo (lat −30.169°, lon −70.806°, 2207 m)
**Noches:** jueves 8 y jueves 15 de octubre de 2026
**Ventana útil:** 21:00 → 01:00 (extendida a 02:00), hora de Chile CLST = UTC−3
**Criterio de visibilidad:** horas sobre 30° de altura (masa de aire < 2) dentro de la ventana, en ambas noches
**Cálculo:** astropy (script `visibilidad_tarea2.py`), contrastado con Staralt
**Campo del MAS500:** 21,6′ × 21,6′ (0,317″/px, medido con Gaia; ver [`verificacion_informe_2025.md`](verificacion_informe_2025.md)). El campo de 24,6′ que se usaba antes venía de una escala errónea

---

## 1. Condiciones de cada noche

| | 8/10 | 15/10 |
|---|---|---|
| Fin crepúsculo astronómico | 21:12 | 21:18 |
| Inicio crepúsculo matutino | 05:52 | 05:42 |
| Iluminación lunar (21:00) | 3 % (casi nueva) | 25 % (creciente) |
| Posición de la Luna | bajo el horizonte toda la noche | AR 17h26m, Dec −27.7° (Sagitario) |
| Altura de la Luna a las 21:00 | −36° | +45° |
| Puesta de la Luna | — | 00:52 |

- **8/10:** la Luna no impone restricciones.
- **15/10:** la Luna está alta al comienzo de la noche y se pone a las 00:52. Se descartaron los candidatos a menos de unos 35° de ella, que en la práctica son los de AR entre 17h y 20h. La banda B es la más afectada por la luz lunar.
- **Staralt** usa "Mean Solar Zone Time" (UTC−4). Su marca de las **24 h corresponde a la 01:00 hora de Chile**.

---

## 2. Objetos de ciencia (10 seleccionados)

Coordenadas J2000 y magnitudes integradas de OpenNGC. "h>30°" son las horas sobre 30° entre las 21:00 y la 01:00. Un valor de 4,0 significa que el objeto está sobre 30° durante toda esa ventana.

| #   | Objeto              | Tipo                  | V   | B    | Tamaño        | h>30° 8/10 | h>30° 15/10 | Alt. máx. | Dist. Luna 15/10 |
| --- | ------------------- | --------------------- | --- | ---- | ------------- | ---------- | ----------- | --------- | ---------------- |
| 1   | NGC 7293 (Helix)    | Nebulosa planetaria   | 7.3 | 7.5  | 16.3′         | 4.0        | 4.0         | 81°       | 67°              |
| 2   | NGC 7009 (Saturno)  | Nebulosa planetaria   | 8.0 | 8.3  | 0.7′          | 4.0        | 4.0         | 71°       | 52°              |
| 3   | NGC 300             | Galaxia Scd           | 8.7 | 8.8  | 19.4′ × 13.1′ | 4.0        | 4.0         | 83°       | 88°              |
| 4   | NGC 7793            | Galaxia Sd            | 9.3 | 9.7  | 10.4′ × 6.0′  | 4.0        | 4.0         | 88°       | 80°              |
| 5   | NGC 247             | Galaxia SABd          | 9.2 | 9.7  | 19.7′ × 5.5′  | 4.0        | 4.0         | 80°       | 95°              |
| 6   | NGC 1097            | Galaxia SBb (Seyfert) | 9.8 | 10.1 | 10.6′ × 6.4′  | 2.4        | 2.9         | 73–79°    | 110°             |
| 7   | NGC 1291            | Galaxia S0/a          | 8.7 | 9.4  | 11.2′ × 9.9′  | 2.2        | 2.6         | 65–70°    | 105°             |
| 8   | NGC 1316 (Fornax A) | Galaxia S0 peculiar   | 8.5 | 9.4  | 13.5′ × 7.7′  | **2.0**    | 2.4         | 65–70°    | 108°             |
| 9   | NGC 1313            | Galaxia SBd           | 9.5 | 9.7  | 11.1′ × 9.1′  | 2.7        | 3.2         | 50–52°    | 83°              |
| 10  | NGC 6744            | Galaxia SABbc         | 9.3 | 9.1  | 15.7′ × 9.8′  | 4.0        | 3.8         | 54°       | **40°**          |

Para los objetos 6 a 9, la altura máxima indicada se alcanza a las 02:00, porque siguen subiendo al cierre de la ventana.

### Encuadre en el campo de 21,6′ × 21,6′

Extensión de cada objeto proyectada en las direcciones este-oeste (E-O) y norte-sur (N-S) del
detector. Se usan los tamaños de la tabla, el ángulo de posición (PA) de SIMBAD y la rotación
medida de la cámara (−0,5°, prácticamente con el norte hacia arriba). "Margen" es el cielo libre
que queda a cada lado en la dirección más justa, con el objeto centrado. "Área" es la fracción del
campo que cubre la elipse del objeto.

| Objeto | PA | Extensión E-O × N-S | Margen mínimo | Área | Encuadre |
|---|---|---|---|---|---|
| NGC 7293 (Helix) | — | 16,3′ × 16,3′ | 2,7′ | 45 % | Cabe; poco cielo libre |
| NGC 7009 | 70° | 0,7′ × 0,7′ | 10,5′ | < 1 % | Sobrado |
| NGC 300 | 114° | 18,5′ × 14,3′ | **1,6′** (E-O) | 43 % | Justo |
| NGC 7793 | 84° | 10,4′ × 6,1′ | 5,6′ | 10 % | Sobrado |
| NGC 247 | 178° | 5,5′ × 19,7′ | **1,0′** (N-S) | 18 % | Muy justo en N-S |
| NGC 1097 | 147° | 7,9′ × 9,5′ | 6,0′ | 11 % | Sobrado |
| NGC 1291 | 156° | 10,1′ × 11,0′ | 5,3′ | 19 % | Sobrado |
| NGC 1316 | 49° | 11,4′ × 10,6′ | 5,1′ | 17 % | Sobrado |
| NGC 1313 | 39° | 9,9′ × 10,4′ | 5,6′ | 17 % | Sobrado |
| NGC 6744 | 16° | 10,4′ × 15,3′ | 3,1′ | 26 % | Cabe |

- Los tamaños corresponden a la isofota de 25 mag/arcsec² en B (D25). Las galaxias tienen luz
  más allá de ese borde, así que el cielo libre real es algo menor que el margen.
- **NGC 247 y NGC 300:** hay que centrarlas bien, porque un error de apuntado de 1′ las deja
  tocando el borde. El fondo se mide en las esquinas y en los márgenes; conviene mencionarlo
  como limitación en el informe.
- **Helix:** el disco principal (16,3′) cabe, pero su halo débil exterior se extiende más allá
  del campo.

### Coordenadas (formato Staralt)

```
NGC7293  22 29 38 -20 50 14
NGC7009  21 04 10 -11 21 47
NGC300   00 54 53 -37 41 03
NGC7793  23 57 49 -32 35 27
NGC247   00 47 08 -20 45 37
NGC1097  02 46 19 -30 16 29
NGC1291  03 17 18 -41 06 29
NGC1316  03 22 41 -37 12 29
NGC1313  03 18 16 -66 29 53
NGC6744  19 09 46 -63 51 27
```

### Observaciones

- **NGC 1316:** está justo en el límite de 2 h el 8/10. Si la ventana se extiende hasta las 02:00, tiene 2,9 h.
- **NGC 6744:** está a 40° de la Luna el 15/10. Conviene observarla al comienzo de la noche, porque además va bajando (alt. 54° a las 21:00).
- **NGC 300, NGC 247 y Helix:** caben en el campo de 21,6′, pero con márgenes de 1,6′, 1,0′ y 2,7′ por lado (ver la tabla de encuadre). Queda poco cielo libre para estimar el fondo, y son objetos de bajo brillo superficial.
- **NGC 1097, 1291, 1316 y 1313:** salen tarde, así que conviene dejarlos para los últimos bloques del plan.

### Reservas

| Objeto         | Tipo                 | V    | B    | Tamaño | h>30° 8/10 | h>30° 15/10 | Dist. Luna 15/10 | Comentario                                                |
| -------------- | -------------------- | ---- | ---- | ------ | ---------- | ----------- | ---------------- | --------------------------------------------------------- |
| NGC 613        | Galaxia SBbc         | 10.4 | 10.7 | 5.5′   | 3.6        | 4.0         | 99°              | Algo débil, bien posicionada                              |
| NGC 1068 (M77) | Galaxia Sb (Seyfert) | 9.3  | 9.7  | 6.1′   | 1.4        | 1.8         | 131°             | Cumple solo si la ventana llega a las 02:00 (2,4 / 2,8 h) |
| NGC 1365       | Galaxia SBb          | 10.1 | 10.4 | 12.0′  | 1.8        | 2.2         | 110°             | Cumple solo si la ventana llega a las 02:00 (2,7 / 3,2 h) |

```
NGC613   01 34 18 -29 25 06
NGC1068  02 42 40 -00 00 47
NGC1365  03 33 36 -36 08 25
```

### Descartados

| Objeto | Motivo |
|---|---|
| NGC 253 | 26.8′ (PA 53°): proyectada mide 21,7′ en E-O, todo el ancho del campo de 21,6′, sin cielo libre |
| NGC 55 | 29.9′ (PA 101°): proyectada mide 29,4′ en E-O, no cabe en el campo de 21,6′ |
| NGC 6822 | 33° de la Luna el 15/10, bajo brillo superficial |
| NGC 6818 | 33° de la Luna el 15/10 |
| NGC 6302 | 10° de la Luna el 15/10, menos de 2 h sobre 30° |
| NGC 628 (M74) | Altura máxima de solo 44° |
| NGC 1360 | 1,5 h sobre 30° antes de la 01:00 el 8/10 |
| NGC 246 | V = 10.9, fuera del rango de magnitud |

---

## 3. Estrellas estándar (8 seleccionadas)

**Criterios:**
- Enanas de tipo A0–A9 V.
- V y B entre 6 y 8, con B = V + (B−V) de Hipparcos.
- Dec entre −22° y −37°, cerca del cenit de Tololo (−30.2°).
- Repartidas en AR para cubrir toda la noche.

**Verificación de variabilidad y multiplicidad:** se revisaron Hipparcos (I/239), WDS (B/wds), VSX (B/vsx) y SIMBAD. Todas las seleccionadas cumplen:
- Aparecen como constantes (HvarType = C) o sin variabilidad detectada en Hipparcos.
- No tienen solución doble (MultFlag vacío) ni entrada en el catálogo de dobles CCDM.
- No aparecen en WDS ni en VSX.
- SIMBAD las clasifica como estrellas normales (sin tipo de variable ni binaria).

| # | Estrella | HIP | AR (J2000) | Dec (J2000) | Tipo | V | B | Culmina 8/10 | Culmina 15/10 | Alt. máx. | Dist. Luna 15/10 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | HD 195500 | 101367 | 20:32:42 | −28:35:44 | A1V | 7.32 | 7.38 | 21:08 | 20:40 | 88° | 39° |
| 2 | HD 202941 | 105318 | 21:19:51 | −27:12:33 | A0V | 7.07 | 7.07 | 21:54 | 21:26 | 87° | 50° |
| 3 | HD 210300 | 109412 | 22:10:00 | −28:17:33 | A5V | 6.44 | 6.59 | 22:44 | 22:16 | 88° | 60° |
| 4 | HD 215863 | 112631 | 22:48:41 | −34:46:23 | A2V | 7.69 | 7.84 | 23:22 | 22:56 | 86° | 66° |
| 5 | HD 220881 | 115796 | 23:27:33 | −27:16:41 | A9V | 7.45 | 7.74 | 00:02 | 23:34 | 87° | 76° |
| 6 | HD 562 | 810 | 00:10:00 | −25:52:29 | A2V | 7.66 | 7.80 | 00:44 | 00:16 | 86° | 86° |
| 7 | HD 8130 | 6257 | 01:20:16 | −36:14:34 | A0V | 7.45 | 7.50 | 01:54 | 01:26 | 84° | 93° |
| 8 | HD 12206 | 9285 | 01:59:20 | −26:25:56 | A0V | 6.79 | 6.81 | 02:34 | 02:06 | 86° | 105° |

Todas pasan más de 3 h sobre 30° dentro de la ventana en ambas noches.

### Coordenadas (formato Staralt)

```
HD195500  20 32 42 -28 35 44
HD202941  21 19 51 -27 12 33
HD210300  22 10 00 -28 17 33
HD215863  22 48 41 -34 46 23
HD220881  23 27 33 -27 16 41
HD562     00 10 00 -25 52 29
HD8130    01 20 16 -36 14 34
HD12206   01 59 20 -26 25 56
```

### Estrellas de reserva (también cumplen todos los criterios)

```
HD207480  21 49 54 -27 24 14
HD212643  22 26 11 -23 40 57
HD223884  23 53 21 -24 13 45
HD225200  00 04 20 -29 16 08
HD7323    01 12 55 -35 44 45
```

### Estrellas descartadas

| Estrella | Motivo |
|---|---|
| HD 222332 (HIP 116750) | Variable en VSX (ASAS J233942-2232.0) |
| HD 6619 (AW Scl) | Variable y binaria espectroscópica |
| HD 7908 (HIP 6108) | Estrella Am (F0VmA3), no es una enana A normal |
| HD 196385, HD 210739, HD 204394 | Compañeras en WDS |
| HD 185404, HD 189830, HD 193281, HD 195206, HD 8487, HD 16087, HD 220455, HD 184439, HD 188113 | Solución doble en Hipparcos / entrada en CCDM |
| HD 182985, HD 186417 | A 26–30° de la Luna el 15/10 |

---

## 4. Verificación con Staralt

Se corrió Staralt (Cerro Tololo, opción "Moon distance") para:
- los 10 objetos el 8/10 y el 15/10;
- las 8 estrellas el 15/10, que es la noche con Luna.

Los resultados coinciden con astropy: iluminación lunar de 2 % y 26 %, Luna en AR 17h43m, Dec −27°51′, y alturas y distancias a la Luna con diferencias de 1–2°.

Para obtener las curvas de visibilidad del informe:
1. Abrir <http://catserver.ing.iac.es/staralt/>.
2. Elegir la fecha, el observatorio "Cerro Tololo Observatory (Chile)" y "Moon distance".
3. Pegar las listas de coordenadas de arriba.
