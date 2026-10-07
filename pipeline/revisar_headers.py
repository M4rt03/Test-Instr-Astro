import glob, os
from astropy.io import fits

claves = ("IMAGETYP", "FILTER", "EXPTIME", "DATE-OBS", "AIRMASS", "EGAIN", "GAIN")
archivos = sorted(glob.glob("datos/2025-10-02/**/*.fit", recursive=True))
print(len(archivos), "archivos encontrados")
for f in archivos[::10]:
    h = fits.getheader(f)
    print(os.path.basename(f), [h.get(k) for k in claves])