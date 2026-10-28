# Pipeline con corrección de transparencia del 15/10 (aparte)

El zeropoint del 15/10 viene de las estándares del 26/10, una noche despejada. Como el 15/10 tuvo
nubes, con ese ZP las magnitudes del 15/10 salen ~0,15 mag más débiles que APASS. Esta carpeta
corre el mismo pipeline con la opción `transparencia`, que mide ese desplazamiento con estrellas
del campo de APASS DR9 y corrige los ZP de la noche. Los resultados quedan en
`con_correccion/resultados`, separados de `pipeline/resultados` (que no cambia).

```
cd pipeline
python run_pipeline.py con_correccion/config.yaml
```

- `config.yaml`: copia de `pipeline/config.yaml` sin la API key (no hace falta en modo `archivo`),
  con rutas a `../datos` y el bloque `transparencia` en la noche del 15/10.
- `apass_ngc7793.ecsv`: 77 estrellas de APASS en 17,3′ alrededor de NGC 7793, consultadas el 16/10
  con `calibrar_campo.py`. Se guardó aquí porque VizieR no respondía el 27/10; si se borra, el
  pipeline vuelve a consultarlo.

## Método

1. Se calibra la noche con el ZP del 26/10 y el término de color de B (c = +2,29).
2. Se cruza `catalogo_BV` con APASS (≤ 2″) y se eligen estrellas no saturadas, con errores < 0,1 mag,
   11 < V < 16 y a más de 6′ de NGC 7793.
3. Δ = mediana con sigma-clipping (2,5σ) de APASS − pipeline, en B y en V por separado (las
   exposiciones de B y V se tomaron en distintos momentos y la nube cambió).
4. ZP_V += Δ_V. En B, para que B final se mueva en Δ_B con el término de color, el ZP de B−V = 0 se
   mueve en Δ_B (1 − c) + c Δ_V.
5. Se vuelven a calcular los catálogos, el diagrama color-magnitud, el mapa de color y el resumen.

## Resultados (`resultados/noche_15oct/transparencia.json`)

| Filtro | Δ (APASS − pipeline) | rms | N | ZP sin corrección | ZP corregido |
|---|---|---|---|---|---|
| V | −0,138 ± 0,021 | 0,065 | 10 de 12 | 21,750 | 21,612 ± 0,021 |
| B | −0,091 ± 0,023 | 0,073 | 10 de 12 | 19,429 (B−V = 0) | 19,229 ± 0,027 |

- La nube quitó 0,14 mag en V y 0,09 mag en B. La diferencia (0,05 ± 0,03) cabe en una nube gris
  que cambió entre las series de B y V.
- El ZP corregido en V (21,61) coincide con el que dieron las estrellas del campo sin estándares
  (21,61–21,63, `extra/calibracion_estrellas_campo`).
- Comparado con APASS (12 estrellas), la mediana pasa de +0,14 a 0,00 en V, de +0,09 a 0,00 en B y
  de −0,05 a 0,00 en B−V. El rms no cambia (0,07 mag): la corrección es un desplazamiento único.
  Dos estrellas débiles (V ≈ 15,4 y 16,0) quedan 0,25–0,43 mag más débiles que APASS y una brillante
  (V = 11,4) 0,15 mag más brillante; el sigma-clipping las deja fuera del Δ.

## Límites

- La mediana queda en cero por construcción: las mismas estrellas sirven para medir el Δ y para
  comprobarlo. No es una validación independiente.
- El Δ también absorbe el error de la extinción supuesta (k_B 0,25, k_V 0,15) y las diferencias de
  escala con APASS (~0,02–0,05 mag). Por eso los ZP corregidos sirven solo para esta noche.
- Con 10 estrellas, el error del Δ (~0,02 mag) domina el error del ZP corregido.
