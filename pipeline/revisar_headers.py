import glob, os
from astropy.io import fits

claves = ("OBJECT", "FILTER", "EXPTIME", "DATE-OBS", "AIRMASS", "CCD-TEMP")
archivos = sorted(glob.glob("datos/2025-10-02/**/*.fit", recursive=True))
print(len(archivos), "archivos encontrados")
for f in archivos:
    h = fits.getheader(f)
    if "HIP" in str(h.get("OBJECT", "")) or "HIP" in f:
        print(os.path.basename(f), [h.get(k) for k in claves])