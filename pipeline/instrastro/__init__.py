"""Pipeline de reducción y análisis para la Tarea 2 de FIS2421 (Instrumentación Astronómica).

Módulos:

- ``io``          : lectura de la configuración y de archivos FITS.
- ``reduction``   : master bias, dark current, master flats y reducción de imágenes.
- ``alignment``   : alineamiento de imágenes por correlación de fase.
- ``astrometry``  : WCS inicial, astrometry.net y resumen de la solución.
- ``photometry``  : estrellas estándar, conversión Tycho -> Johnson, zeropoints.
- ``detection``   : fondo, detección, catálogos con dos configuraciones de deblending.
- ``seeing``      : FWHM de estrellas por ajuste gaussiano, límite de difracción, variación de la PSF.
- ``color``       : imagen RGB y mapa de color B-V.
- ``planning``    : curvas de visibilidad y bitácora de observaciones.
- ``pipeline``    : orquesta todos los pasos para cada noche.
- ``synthetic``   : genera un set de datos sintético para probar el pipeline.
"""

__version__ = "1.0.0"
