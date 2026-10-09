# Test-Instr-Astro

Material para la Tarea 2 de FIS2421 Instrumentación Astronómica (observaciones en B y V con el MAS500).

- `pipeline/`: pipeline limpio y probado para procesar tus propios datos: reducción, astrometría, fotometría, catálogos, seeing e imagen a color. Instrucciones en [`pipeline/README.md`](pipeline/README.md).
- `planificacion/candidatos_2026.md`: objetos candidatos y estrellas estándar para las noches del 8 y 15 de octubre de 2026, con visibilidad, distancia a la Luna y verificación en Staralt.
- `planificacion/tiempos_exposicion.md`: cálculo de tiempos de exposición (ecuación de S/N, saturación, plan por bloque) y la calculadora `tiempos_exposicion.py`; el efecto de la Luna en las noches del 15 y 22/10 se calcula con `cielo_luna.py`.
- `planificacion/calibraciones.md`: plan y tiempos de bias, darks, flats y estrellas estándar; las magnitudes Johnson B y V de las estándares se consultan con astroquery en `planificacion/estrellas_estandar.py`.
- `planificacion/pendientes.md`: lista de lo que hay que confirmar antes, durante y después de las observaciones.
- `planificacion/verificacion_informe_2025.md`: datos de la cámara según los headers, escala de placa medida con Gaia (`medir_escala.py`) y revisión de los cálculos del informe de 2025.
- Raíz del repositorio: material de referencia de 2025 (enunciado, informe sobre NGC 6618, notebooks originales y apuntes).
