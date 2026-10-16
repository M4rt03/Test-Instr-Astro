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

Las tablas quedan en `resultados/`. El ZP sale a la masa de aire del apilado; para pasarlo a
fuera de la atmósfera se usa un k supuesto (B 0,25, V 0,15), porque sin estándares a dos masas
de aire no se puede medir.
