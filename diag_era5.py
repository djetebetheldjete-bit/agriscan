# diag_era5.py — Diagnostic ERA5 : verifie ce que le dataset repond vraiment
import ee

PROJET = "ambient-fuze-428510-q0"
CLE = "cle-service.json"
import os
if os.path.exists(CLE):
    ee.Initialize(ee.ServiceAccountCredentials(None, CLE), project=PROJET)
else:
    ee.Initialize(project=PROJET)
print("GEE OK\n")

# Zone de test (Sinfra)
zone = ee.Geometry.Rectangle([-5.945, 6.835, -5.895, 6.875])

print("=== TEST 1 : le dataset ERA5_LAND/DAILY_AGGR existe-t-il et a-t-il des donnees 2024 ? ===")
for periode in [("2024-11-01","2024-11-30"), ("2024-12-01","2024-12-15"),
                ("2024-12-15","2024-12-31"), ("2024-06-01","2024-06-30")]:
    col = (ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
           .filterDate(periode[0], periode[1]).filterBounds(zone))
    try:
        n = col.size().getInfo()
        print(f"  {periode[0]} -> {periode[1]} : {n} images")
    except Exception as e:
        print(f"  {periode[0]} -> {periode[1]} : ERREUR {e}")

print("\n=== TEST 2 : quelles bandes contient une image ERA5 ? ===")
try:
    col = (ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
           .filterDate("2024-06-01","2024-06-30").filterBounds(zone))
    if col.size().getInfo() > 0:
        bandes = col.first().bandNames().getInfo()
        print("  Bandes disponibles (extrait) :")
        # afficher celles qu'on utilise
        cherchees = ["temperature_2m_max","temperature_2m_min","dewpoint_temperature_2m",
                     "u_component_of_wind_10m","v_component_of_wind_10m","surface_pressure",
                     "surface_net_solar_radiation_sum","surface_net_thermal_radiation_sum",
                     "temperature_2m"]
        for b in cherchees:
            present = "OK" if b in bandes else "ABSENTE !!!"
            print(f"    {b} : {present}")
    else:
        print("  Pas d'image pour extraire les bandes")
except Exception as e:
    print(f"  ERREUR : {e}")

print("\n=== TEST 3 : calcul ET0 reel sur juin 2024 (valeur attendue ~3-5 mm/j) ===")
try:
    import math as _math
    col = (ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
           .filterDate("2024-06-01","2024-06-30").filterBounds(zone))
    def et0_image(img):
        tmax = img.select("temperature_2m_max").subtract(273.15)
        tmin = img.select("temperature_2m_min").subtract(273.15)
        tmean = tmax.add(tmin).divide(2)
        tdew = img.select("dewpoint_temperature_2m").subtract(273.15)
        u10 = img.select("u_component_of_wind_10m")
        v10 = img.select("v_component_of_wind_10m")
        vent10 = u10.hypot(v10)
        u2 = vent10.multiply(4.87).divide(_math.log(67.8*10-5.42))
        def e0(t): return t.multiply(17.27).divide(t.add(237.3)).exp().multiply(0.6108)
        es = e0(tmax).add(e0(tmin)).divide(2)
        ea = e0(tdew)
        delta = e0(tmean).multiply(4098).divide(tmean.add(237.3).pow(2))
        P = img.select("surface_pressure").divide(1000)
        gamma = P.multiply(0.000665)
        rns = img.select("surface_net_solar_radiation_sum").divide(1e6)
        rnl = img.select("surface_net_thermal_radiation_sum").divide(1e6)
        Rn = rns.add(rnl)
        num = (delta.multiply(0.408).multiply(Rn)
               .add(gamma.multiply(tmean.add(273).pow(-1).multiply(900)).multiply(u2).multiply(es.subtract(ea))))
        den = delta.add(gamma.multiply(u2.multiply(0.34).add(1)))
        return num.divide(den).rename("ET0")
    et0 = col.map(et0_image).mean().clip(zone)
    val = et0.reduceRegion(reducer=ee.Reducer.mean(), geometry=zone, scale=1000, maxPixels=1e9, bestEffort=True).get("ET0").getInfo()
    print(f"  ET0 moyen juin 2024 = {val}")
except Exception as e:
    print(f"  ERREUR calcul ET0 : {e}")


print("\n=== TEST 4 : reproduction EXACTE de la route /irrigation ===")
import math as _math
from datetime import datetime, timedelta

def calcul_et0_gee(zone, date_debut, date_fin):
    col = (ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
           .filterDate(date_debut, date_fin).filterBounds(zone))
    def et0_image(img):
        tmax = img.select("temperature_2m_max").subtract(273.15)
        tmin = img.select("temperature_2m_min").subtract(273.15)
        tmean = tmax.add(tmin).divide(2)
        tdew = img.select("dewpoint_temperature_2m").subtract(273.15)
        u10 = img.select("u_component_of_wind_10m"); v10 = img.select("v_component_of_wind_10m")
        vent10 = u10.hypot(v10)
        u2 = vent10.multiply(4.87).divide(_math.log(67.8*10-5.42))
        def e0(t): return t.multiply(17.27).divide(t.add(237.3)).exp().multiply(0.6108)
        es = e0(tmax).add(e0(tmin)).divide(2); ea = e0(tdew)
        delta = e0(tmean).multiply(4098).divide(tmean.add(237.3).pow(2))
        P = img.select("surface_pressure").divide(1000); gamma = P.multiply(0.000665)
        rns = img.select("surface_net_solar_radiation_sum").divide(1e6)
        rnl = img.select("surface_net_thermal_radiation_sum").divide(1e6)
        Rn = rns.add(rnl)
        num = (delta.multiply(0.408).multiply(Rn)
               .add(gamma.multiply(tmean.add(273).pow(-1).multiply(900)).multiply(u2).multiply(es.subtract(ea))))
        den = delta.add(gamma.multiply(u2.multiply(0.34).add(1)))
        return num.divide(den).rename("ET0").copyProperties(img, ["system:time_start"])
    return col.map(et0_image).mean().clip(zone)

def lire_et0_moyen(et0_img, zone):
    for sc in (1000, 5000, 11132):
        val = et0_img.reduceRegion(reducer=ee.Reducer.mean(), geometry=zone,
                                   scale=sc, maxPixels=1e9, bestEffort=True).get("ET0").getInfo()
        if val is not None:
            return val
    centre = zone.centroid(maxError=1)
    return et0_img.reduceRegion(reducer=ee.Reducer.first(), geometry=centre,
                                scale=11132, maxPixels=1e9, bestEffort=True).get("ET0").getInfo()

date_fin = "2024-12-15"
fin_dt = datetime.strptime(date_fin, "%Y-%m-%d")
debut_dt = fin_dt - timedelta(days=60)
print(f"  Fenetre : {debut_dt.strftime('%Y-%m-%d')} -> {date_fin}")
col_test = (ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
            .filterDate(debut_dt.strftime("%Y-%m-%d"), date_fin).filterBounds(zone))
print(f"  n_images = {col_test.size().getInfo()}")
et0_img = calcul_et0_gee(zone, debut_dt.strftime("%Y-%m-%d"), date_fin)
val = lire_et0_moyen(et0_img, zone)
print(f"  ET0 via lire_et0_moyen = {val}")
if val is None:
    print("  >>> REPRODUIT LE BUG : lire_et0_moyen renvoie None")
else:
    print("  >>> PAS DE BUG ICI : le calcul marche. Le probleme est ailleurs (passerelle pas rechargee).")
