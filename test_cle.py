# test_cle.py — Verifie que le compte de service GEE fonctionne
# Objectif : prouver que le "robot" (compte de service) peut initialiser
# Earth Engine avec sa cle JSON, sans authentification personnelle.

import ee

# Chemin de la cle du compte de service (dans le meme dossier)
CLE = "cle-service.json"
PROJET = "ambient-fuze-428510-q0"

print("1) Lecture de la cle et connexion a Earth Engine...")
try:
    credentials = ee.ServiceAccountCredentials(None, CLE)
    ee.Initialize(credentials, project=PROJET)
    print("   -> Earth Engine initialise avec le compte de service. OK")
except Exception as e:
    print("   -> ECHEC de l'initialisation :")
    print("   ", e)
    raise SystemExit(1)

print("2) Petit calcul de test (NDVI moyen sur une mini-zone Sinfra)...")
try:
    zone = ee.Geometry.Rectangle([-5.945, 6.835, -5.895, 6.875])
    col = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
           .filterBounds(zone)
           .filterDate("2024-01-01", "2024-12-31")
           .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 20)))
    ndvi = col.median().normalizedDifference(["B8", "B4"])
    moyenne = ndvi.reduceRegion(
        reducer=ee.Reducer.mean(), geometry=zone, scale=10, maxPixels=1e9
    ).get("nd").getInfo()
    print("   -> NDVI moyen calcule :", round(moyenne, 3) if moyenne is not None else None)
    print()
    print("SUCCES : le compte de service fonctionne. On peut deployer en ligne.")
except Exception as e:
    print("   -> Le calcul a echoue :")
    print("   ", e)
    raise SystemExit(1)
