# Datos del 15/10/2026 con el método de 2025 (fuera de la tarea)

`reduccion_2025_15oct.py` reduce NGC 7793 (noche del 15/10/2026) siguiendo los pasos de los
notebooks del grupo de 2025 (`Reducción3.333.ipynb` y `reduccion y det. de objetos.ipynb`), para
ver cómo se ven los datos del jueves procesados así y compararlos con el pipeline actual.

```
cd extra\pipeline_2025
python reduccion_2025_15oct.py
```

## Qué replica y qué hubo que cambiar

| Paso de 2025 | Aquí |
|---|---|
| Master bias = promedio de 10 bias | No hay bias en 2026: se usan los 5 darks de 120 s como nivel cero |
| Corriente oscura = promedio(dark − bias)/100 s | Igual, con los 5 darks de 240 s |
| Master dark de 30 s = bias + corriente × 30 | Igual, para 120 s |
| Flats: flat − bias − corriente × t, promedio, normalizado por el promedio | Igual |
| Ciencia: (cruda − master dark)/flat y promedio, sin alinear ni excluir | Igual: 5 imágenes en B y 4 en V de 120 s, incluidas las 2 con estrellas dobles |
| × BPM.fit | El BPM de 2025 no está en el repo: se omite |
| PIXSCALE = 0,36 y WCS desde OBJCTRA/OBJCTDEC | Igual |
| Detección: nsigma = 1, npixels = 10, deblend con contrast = 1 | Igual |

Las magnitudes de 2025 (ZP inventado de ~39 y `np.log` en vez de `log10`) no se replican.

## Resultado

- **`reducidas_metodo2025_15oct.png`:** las imágenes B y V reducidas, con la misma escala zscale
  del notebook. La galaxia se ve, pero las esquinas quedan ruidosas: al normalizar el flat por el
  promedio (no por la mediana del centro) y sin enmascarar las esquinas viñeteadas, el flat divide
  por valores muy bajos ahí.
- **`comparacion_estrella_15oct.png`:** recorte de 1,3′ alrededor de una estrella de V ≈ 11,4. Con
  el método de 2025 la estrella sale **doble en B y triple/alargada en V**: las imágenes se
  promedian sin alinear y el telescopio se movió hasta 50 px entre exposiciones (además de las dos
  imágenes con estrellas dobles por el viento). En el pipeline actual, que alinea con las estrellas
  y excluye esas imágenes, sale una sola.
- **Detección:** con nsigma = 1 se detectan ~10 000 (V) y ~12 500 (B) "fuentes", casi todas ruido, y
  la galaxia queda unida al fondo en un solo segmento de ~600 000 px. El catálogo del pipeline actual
  (umbral de 3σ, fondo local) tiene ~1300 fuentes por filtro.
- **Astrometría:** el WCS de 2025 queda centrado a 40′ de NGC 7793 (el `OBJCTRA` de esa noche estaba
  corrido) y con escala 0,36″/px en vez de 0,316″/px.

Los FITS y catálogos generados quedan en `salida/` (no se suben a git).
