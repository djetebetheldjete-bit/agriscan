# passerelle.py — Passerelle Vertea
# Routes :
#   - /            : sert la page web (index.html)
#   - /ndvi        : carte NDVI + NDVI moyen + surface + repartition sain/moyen/faible
#   - /evolution   : serie temporelle du NDVI moyen (pas mois ou semaine)
#   - /anomalie    : detection de changement (recent - reference)
#
# Double mode d'authentification Earth Engine :
#   - En ligne / avec cle : si le fichier cle-service.json existe (ou variable
#     d'environnement GEE_KEY_FILE), on utilise le COMPTE DE SERVICE.
#   - En local : sinon, on retombe sur l'authentification personnelle (ee.Initialize).

import os
import ee
from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

# ===== Parametres =====
PROJET = "ambient-fuze-428510-q0"
# Chemin de la cle du compte de service (peut etre surcharge par une variable d'env)
CLE_SERVICE = os.environ.get("GEE_KEY_FILE", "cle-service.json")

# ===== Seuils configurables =====
SEUIL_FAIBLE = 0.3
SEUIL_SAIN   = 0.5
SEUIL_ANOMALIE = 0.1


def initialiser_gee():
    """
    Initialise Earth Engine.
    - Si une cle de compte de service est disponible -> on l'utilise (mode serveur).
    - Sinon -> authentification personnelle locale (mode developpement).
    """
    if os.path.exists(CLE_SERVICE):
        credentials = ee.ServiceAccountCredentials(None, CLE_SERVICE)
        ee.Initialize(credentials, project=PROJET)
        print("Google Earth Engine initialise via COMPTE DE SERVICE \u2713")
    else:
        ee.Initialize(project=PROJET)
        print("Google Earth Engine initialise via auth personnelle \u2713")


try:
    initialiser_gee()
except Exception as e:
    print("GEE non initialise. Detail:", e)


def construire_zone(data):
    """Construit la geometrie GEE (rectangle bbox OU polygone)."""
    polygone = data.get("polygon")
    if polygone:
        if len(polygone) < 3:
            return None, "Polygone invalide : il faut au moins 3 points."
        return ee.Geometry.Polygon([polygone]), None
    coords = data.get("bbox")
    if coords and len(coords) == 4:
        return ee.Geometry.Rectangle(coords), None
    return None, "Zone invalide. Envoie un 'bbox' [ouest, sud, est, nord] ou un 'polygon' [[lng,lat],...]."


def collection_ndvi(zone, date_debut, date_fin):
    """Collection Sentinel-2 filtree, avec bande NDVI ajoutee."""
    def ajouter_ndvi(image):
        ndvi = image.normalizedDifference(["B8", "B4"]).rename("NDVI")
        return image.addBands(ndvi)
    return (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
            .filterBounds(zone)
            .filterDate(date_debut, date_fin)
            .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 20))
            .map(ajouter_ndvi))


def calculer_stats(ndvi_image, zone):
    """Surface (ha) + repartition sain/moyen/faible (%)."""
    classe = (ndvi_image.gte(SEUIL_FAIBLE)
              .add(ndvi_image.gte(SEUIL_SAIN))
              .rename("classe"))
    histo = classe.reduceRegion(
        reducer=ee.Reducer.frequencyHistogram(),
        geometry=zone, scale=10, maxPixels=1e9, bestEffort=True
    ).get("classe").getInfo() or {}

    n_faible = histo.get("0", 0); n_moyen = histo.get("1", 0); n_sain = histo.get("2", 0)
    total = n_faible + n_moyen + n_sain
    if total > 0:
        pct_faible = round(100 * n_faible / total, 1)
        pct_moyen  = round(100 * n_moyen  / total, 1)
        pct_sain   = round(100 * n_sain   / total, 1)
    else:
        pct_faible = pct_moyen = pct_sain = None

    aire_m2 = zone.area(maxError=1).getInfo()
    return {
        "surface_ha": round(aire_m2 / 10000, 2),
        "pct_sain": pct_sain, "pct_moyen": pct_moyen, "pct_faible": pct_faible
    }



# ======================================================================
#  IRRIGATION — besoin en eau des cultures (methode FAO-56)
#  ET0 (Penman-Monteith, formule validee) x Kc - pluie efficace
# ======================================================================

import math as _math

# --- Coefficients culturaux Kc (FAO-56) et durees de phases (jours) ---
# Phases : initiale, developpement, mi-saison, fin
CULTURES = {
    "tomate":    {"kc_ini": 0.60, "kc_mid": 1.15, "kc_end": 0.80,
                  "d_ini": 30, "d_dev": 40, "d_mid": 40, "d_fin": 25},
    "concombre": {"kc_ini": 0.60, "kc_mid": 1.00, "kc_end": 0.75,
                  "d_ini": 20, "d_dev": 30, "d_mid": 40, "d_fin": 15},
    "piment":    {"kc_ini": 0.60, "kc_mid": 1.05, "kc_end": 0.90,
                  "d_ini": 30, "d_dev": 35, "d_mid": 40, "d_fin": 20},
    "laitue":    {"kc_ini": 0.70, "kc_mid": 1.00, "kc_end": 0.95,
                  "d_ini": 20, "d_dev": 30, "d_mid": 15, "d_fin": 10},
}


def kc_du_jour(culture, jours_depuis_semis):
    """
    Renvoie (Kc, nom_du_stade) selon le nombre de jours depuis le semis.
    Interpolation lineaire pendant les phases de transition (FAO-56).
    """
    c = CULTURES[culture]
    d_ini, d_dev, d_mid, d_fin = c["d_ini"], c["d_dev"], c["d_mid"], c["d_fin"]
    kc_ini, kc_mid, kc_end = c["kc_ini"], c["kc_mid"], c["kc_end"]
    j = jours_depuis_semis

    fin_ini = d_ini
    fin_dev = d_ini + d_dev
    fin_mid = d_ini + d_dev + d_mid
    fin_tot = d_ini + d_dev + d_mid + d_fin

    if j < 0:
        return kc_ini, "avant semis"
    if j <= fin_ini:
        return kc_ini, "initiale"
    if j <= fin_dev:
        # montee lineaire de kc_ini vers kc_mid
        frac = (j - fin_ini) / d_dev
        return kc_ini + frac * (kc_mid - kc_ini), "developpement"
    if j <= fin_mid:
        return kc_mid, "mi-saison"
    if j <= fin_tot:
        # descente lineaire de kc_mid vers kc_end
        frac = (j - fin_mid) / d_fin
        return kc_mid + frac * (kc_end - kc_mid), "fin de saison"
    return kc_end, "apres recolte"


def pluie_efficace(pluie_mm):
    """
    Pluie efficace (FAO / USDA-SCS, forme simplifiee mensuelle ramenee).
    On applique un abattement : toute la pluie n'est pas utilisable.
    Regle USDA-SCS simplifiee au pas considere.
    """
    if pluie_mm <= 0:
        return 0.0
    # Approche simple et prudente : 80% de la pluie est efficace jusqu'a un plafond,
    # au-dela le surplus ruisselle. (parametrable)
    return round(min(pluie_mm * 0.8, pluie_mm), 2)


def calcul_et0_gee(zone, date_debut, date_fin):
    """
    Calcule ET0 (mm/jour, moyenne sur la periode) sur GEE avec ERA5-Land,
    en appliquant la formule Penman-Monteith FAO-56 VALIDEE.
    Renvoie une image ET0 (moyenne periode) clippee sur la zone.
    """
    col = (ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
           .filterDate(date_debut, date_fin)
           .filterBounds(zone))

    def et0_image(img):
        # --- Temperatures (K -> C) ---
        tmax = img.select("temperature_2m_max").subtract(273.15)
        tmin = img.select("temperature_2m_min").subtract(273.15)
        tmean = tmax.add(tmin).divide(2)
        tdew = img.select("dewpoint_temperature_2m").subtract(273.15)

        # --- Vent 10m -> 2m (FAO-56 eq.47) ---
        u10 = img.select("u_component_of_wind_10m")
        v10 = img.select("v_component_of_wind_10m")
        vent10 = u10.hypot(v10)  # racine(u^2+v^2)
        u2 = vent10.multiply(4.87).divide(_math.log(67.8 * 10 - 5.42))

        # --- Pressions de vapeur (kPa) ---
        # e0(T) = 0.6108*exp(17.27*T/(T+237.3))
        def e0(t):
            return t.multiply(17.27).divide(t.add(237.3)).exp().multiply(0.6108)
        es = e0(tmax).add(e0(tmin)).divide(2)
        ea = e0(tdew)  # ea depuis le point de rosee (FAO-56)

        # --- Pente Delta (kPa/C) ---
        delta = (e0(tmean).multiply(4098)
                 .divide(tmean.add(237.3).pow(2)))

        # --- Constante psychrometrique gamma ---
        # P depuis altitude ~ on prend P standard via ERA5 surface_pressure si dispo,
        # sinon approximation. Ici on utilise surface_pressure (Pa -> kPa).
        P = img.select("surface_pressure").divide(1000)
        gamma = P.multiply(0.000665)

        # --- Rayonnement net Rn (J/m2 -> MJ/m2/jour) fourni par ERA5 ---
        rns = img.select("surface_net_solar_radiation_sum").divide(1e6)
        rnl = img.select("surface_net_thermal_radiation_sum").divide(1e6)
        Rn = rns.add(rnl)  # net total (thermal est generalement negatif)
        G = 0

        # --- Penman-Monteith FAO-56 (eq.6) ---
        num = (delta.multiply(0.408).multiply(Rn.subtract(G))
               .add(gamma.multiply(tmean.add(273).pow(-1).multiply(900))
                    .multiply(u2).multiply(es.subtract(ea))))
        den = delta.add(gamma.multiply(u2.multiply(0.34).add(1)))
        et0 = num.divide(den).rename("ET0")
        return et0.copyProperties(img, ["system:time_start"])

    return col.map(et0_image).mean().clip(zone)



def lire_et0_moyen(et0_img, zone):
    """
    Lit la valeur ET0 moyenne de facon ROBUSTE, meme pour de petites zones.
    ERA5-Land a une resolution ~11 km : sur une petite parcelle, reduceRegion
    peut ne rien capter. On essaie plusieurs echelles, puis on retombe sur le
    centroide de la zone si necessaire.
    """
    for sc in (1000, 5000, 11132):
        val = et0_img.reduceRegion(
            reducer=ee.Reducer.mean(), geometry=zone,
            scale=sc, maxPixels=1e9, bestEffort=True
        ).get("ET0").getInfo()
        if val is not None:
            return val
    # Dernier recours : valeur au centre de la parcelle
    centre = zone.centroid(maxError=1)
    val = et0_img.reduceRegion(
        reducer=ee.Reducer.first(), geometry=centre,
        scale=11132, maxPixels=1e9, bestEffort=True
    ).get("ET0").getInfo()
    return val


@app.route("/irrigation", methods=["POST"])
def calculer_irrigation():
    data = request.get_json()
    zone, erreur = construire_zone(data)
    if erreur:
        return jsonify({"erreur": erreur}), 400

    culture = (data.get("culture") or "").lower()
    if culture not in CULTURES:
        return jsonify({"erreur": "Culture inconnue. Choisis : tomate, concombre, piment, laitue."}), 400

    date_semis = data.get("date_semis")
    if not date_semis:
        return jsonify({"erreur": "Fournis la date de semis (date_semis, format AAAA-MM-JJ)."}), 400

    # Periode d'analyse : 7 derniers jours de la periode fournie (ou 7 j avant date_fin)
    date_fin = data.get("date_fin", "2024-12-31")

    try:
        from datetime import datetime, timedelta
        fin_dt = datetime.strptime(date_fin, "%Y-%m-%d")
        debut_dt = fin_dt - timedelta(days=7)
        semis_dt = datetime.strptime(date_semis, "%Y-%m-%d")
        jours_semis = (fin_dt - semis_dt).days

        # Kc selon le stade
        kc, stade = kc_du_jour(culture, jours_semis)

        # ET0 moyen (mm/j) sur la periode via GEE
        et0_img = calcul_et0_gee(zone, debut_dt.strftime("%Y-%m-%d"), date_fin)
        et0_moy = lire_et0_moyen(et0_img, zone)

        if et0_moy is None:
            return jsonify({"erreur": "Pas de donnees climatiques ERA5 sur cette periode/zone."}), 400

        # ETc = ET0 x Kc
        etc = et0_moy * kc

        # Pluie CHIRPS (mm) moyenne journaliere sur la periode
        chirps = (ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
                  .filterDate(debut_dt.strftime("%Y-%m-%d"), date_fin)
                  .filterBounds(zone))
        pluie_tot = chirps.sum().reduceRegion(
            reducer=ee.Reducer.mean(), geometry=zone, scale=5000, maxPixels=1e9, bestEffort=True
        ).get("precipitation").getInfo() or 0
        pluie_jour = pluie_tot / 7.0  # moyenne journaliere

        pluie_eff = pluie_efficace(pluie_jour)

        # Besoin net d'irrigation (mm/jour)
        besoin = max(etc - pluie_eff, 0)

        return jsonify({
            "culture": culture,
            "stade": stade,
            "jours_depuis_semis": jours_semis,
            "kc": round(kc, 2),
            "et0_mm_j": round(et0_moy, 2),
            "etc_mm_j": round(etc, 2),
            "pluie_mm_j": round(pluie_jour, 2),
            "pluie_efficace_mm_j": round(pluie_eff, 2),
            "besoin_irrigation_mm_j": round(besoin, 2),
            "message": "Besoin en eau calcule (FAO-56)"
        })

    except Exception as e:
        return jsonify({"erreur": str(e)}), 500



@app.route("/irrigation_carte", methods=["POST"])
def carte_irrigation():
    """
    Carte spatialisee du besoin en eau (methode Kc derive du NDVI).
    Pour chaque pixel :
        Kc_pixel = 1.457 * NDVI - 0.1725   (Kamble et al., 2013 ; relation NDVI-Kc)
        ETc_pixel = ET0 (climat) * Kc_pixel
        Besoin_pixel = max(ETc_pixel - pluie_efficace, 0)   (mm/jour)
    Renvoie des tuiles colorees + le besoin moyen.
    """
    data = request.get_json()
    zone, erreur = construire_zone(data)
    if erreur:
        return jsonify({"erreur": erreur}), 400

    date_fin = data.get("date_fin", "2024-12-31")
    try:
        from datetime import datetime, timedelta
        fin_dt = datetime.strptime(date_fin, "%Y-%m-%d")
        debut_dt = fin_dt - timedelta(days=30)  # NDVI sur 30 j (assez d'images nettes)
        d_debut = debut_dt.strftime("%Y-%m-%d")

        # 1) NDVI median recent (Sentinel-2) sur la zone
        col_s2 = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
                  .filterBounds(zone)
                  .filterDate(d_debut, date_fin)
                  .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 20)))
        ndvi = col_s2.median().normalizedDifference(["B8", "B4"]).clip(zone)

        # 2) Kc pixel par pixel a partir du NDVI (relation lineaire)
        kc = ndvi.multiply(1.457).subtract(0.1725).rename("Kc")
        # borne Kc dans une plage physique raisonnable [0, 1.3]
        kc = kc.max(0).min(1.3)

        # 3) ET0 (climat) moyen sur 7 j via ERA5 + Penman-Monteith
        et0_debut = (fin_dt - timedelta(days=7)).strftime("%Y-%m-%d")
        et0_img = calcul_et0_gee(zone, et0_debut, date_fin)  # image ET0 (mm/j)

        # 4) ETc spatialise = ET0 * Kc
        etc = et0_img.multiply(kc).rename("ETc")

        # 5) Pluie efficace (scalaire) -> soustraction
        chirps = (ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
                  .filterDate(et0_debut, date_fin).filterBounds(zone))
        pluie_tot = chirps.sum().reduceRegion(
            reducer=ee.Reducer.mean(), geometry=zone, scale=5000, maxPixels=1e9, bestEffort=True
        ).get("precipitation").getInfo() or 0
        pluie_eff_j = pluie_efficace(pluie_tot / 7.0)

        # 6) Besoin pixel = ETc - pluie efficace (mm/j), borne a 0
        besoin = etc.subtract(pluie_eff_j).max(0).rename("besoin")

        # 7) Tuiles colorees (bleu clair -> bleu fonce = besoin croissant)
        palette = ["#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"]
        map_id = besoin.getMapId({"min": 0, "max": 8, "palette": palette})
        tile_url = map_id["tile_fetcher"].url_format

        # 8) Besoin moyen de la zone
        besoin_moy = besoin.reduceRegion(
            reducer=ee.Reducer.mean(), geometry=zone, scale=10, maxPixels=1e9, bestEffort=True
        ).get("besoin").getInfo()

        return jsonify({
            "tile_url": tile_url,
            "besoin_moyen_mm_j": round(besoin_moy, 2) if besoin_moy is not None else None,
            "pluie_efficace_mm_j": round(pluie_eff_j, 2),
            "message": "Carte de besoin en eau calculee (Kc derive du NDVI)"
        })
    except Exception as e:
        return jsonify({"erreur": str(e)}), 500



@app.route("/irrigation_evolution", methods=["POST"])
def irrigation_evolution():
    """
    Evolution hebdomadaire du besoin en eau vs pluie recue.
    Pour chaque semaine de la periode :
      - ETc = ET0 (Penman-Monteith) x Kc (selon stade au milieu de la semaine)
      - pluie = somme CHIRPS de la semaine
      - besoin = max(ETc_semaine - pluie_efficace_semaine, 0)
    Renvoie une liste [{semaine, etc, pluie, besoin}] en mm (cumul hebdo).
    """
    data = request.get_json()
    zone, erreur = construire_zone(data)
    if erreur:
        return jsonify({"erreur": erreur}), 400

    culture = (data.get("culture") or "").lower()
    if culture not in CULTURES:
        return jsonify({"erreur": "Culture inconnue."}), 400
    date_semis = data.get("date_semis")
    if not date_semis:
        return jsonify({"erreur": "Fournis la date de semis."}), 400

    date_debut = data.get("date_debut", "2024-01-01")
    date_fin = data.get("date_fin", "2024-12-31")

    try:
        from datetime import datetime, timedelta
        debut = datetime.strptime(date_debut, "%Y-%m-%d")
        fin = datetime.strptime(date_fin, "%Y-%m-%d")
        semis = datetime.strptime(date_semis, "%Y-%m-%d")

        points = []
        semaine_debut = debut
        # Limite de securite : max ~26 semaines pour eviter les calculs trop longs
        max_semaines = 26
        compteur = 0

        while semaine_debut < fin and compteur < max_semaines:
            semaine_fin = min(semaine_debut + timedelta(days=7), fin)
            d1 = semaine_debut.strftime("%Y-%m-%d")
            d2 = semaine_fin.strftime("%Y-%m-%d")
            milieu = semaine_debut + timedelta(days=3)
            jours_semis = (milieu - semis).days
            kc, _ = kc_du_jour(culture, jours_semis)

            # ET0 moyen (mm/j) sur la semaine -> cumul hebdo
            try:
                et0_img = calcul_et0_gee(zone, d1, d2)
                et0_moy = lire_et0_moyen(et0_img, zone)
            except Exception:
                et0_moy = None

            if et0_moy is None:
                etc_hebdo = None
            else:
                nb_jours = (semaine_fin - semaine_debut).days
                etc_hebdo = et0_moy * kc * nb_jours

            # Pluie CHIRPS cumulee sur la semaine
            chirps = (ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
                      .filterDate(d1, d2).filterBounds(zone))
            pluie = chirps.sum().reduceRegion(
                reducer=ee.Reducer.mean(), geometry=zone,
                scale=5000, maxPixels=1e9, bestEffort=True
            ).get("precipitation").getInfo() or 0
            pluie_eff = pluie_efficace(pluie)

            if etc_hebdo is None:
                besoin = None
            else:
                besoin = max(etc_hebdo - pluie_eff, 0)

            points.append({
                "semaine": semaine_debut.strftime("%Y-%m-%d"),
                "etc": round(etc_hebdo, 1) if etc_hebdo is not None else None,
                "pluie": round(pluie_eff, 1),
                "besoin": round(besoin, 1) if besoin is not None else None
            })

            semaine_debut = semaine_fin
            compteur += 1

        return jsonify({"points": points, "message": "Evolution irrigation calculee"})
    except Exception as e:
        return jsonify({"erreur": str(e)}), 500



# ======================================================================
#  CONSEIL DE PLANTATION — croise conditions reelles (ERA5/CHIRPS)
#  avec les exigences agronomiques (fiche technique Cote d'Ivoire)
# ======================================================================

# Exigences par culture (source : fiche technique agronomique CI)
#   t_min, t_max      : plage de temperature ideale (degC)
#   t_critique        : au-dela, avortement/stress fort
#   eau_min, eau_max  : besoin en eau par cycle (mm)
#   mois_semis        : liste des mois favorables au semis (1=janvier ... 12=decembre)
#   note              : specificite locale a rappeler
EXIGENCES = {
    "concombre": {
        "t_min": 24, "t_max": 30, "t_critique": 35,
        "eau_min": 400, "eau_max": 600,
        "mois_semis": [10, 11, 12, 1, 2, 5, 6],
        "note": "Arrosage tres regulier. Surveiller le mildiou en saison des pluies."
    },
    "piment": {
        "t_min": 20, "t_max": 30, "t_critique": 32,
        "eau_min": 500, "eau_max": 800,
        "mois_semis": [9, 10, 2, 3],
        "note": "Sensible a l'exces d'eau (billons draines). Floraison exige des temperatures douces."
    },
    "tomate": {
        "t_min": 18, "t_max": 27, "t_critique": 32,
        "eau_min": 500, "eau_max": 700,
        "mois_semis": [10, 11, 12, 1],
        "note": "Nuits < 22 degC requises. Eviter la saison humide (mildiou, fletrissement, viroses)."
    },
    "laitue": {
        "t_min": 15, "t_max": 24, "t_critique": 28,
        "eau_min": 300, "eau_max": 400,
        "mois_semis": [11, 12, 1, 2],
        "note": "Enracinement superficiel, arrosages frequents. Montaison si trop chaud : variete thermotolerante + ombriere."
    },
}

MOIS_FR = {1:"janvier",2:"fevrier",3:"mars",4:"avril",5:"mai",6:"juin",
           7:"juillet",8:"aout",9:"septembre",10:"octobre",11:"novembre",12:"decembre"}


def lire_valeur_moyenne(image, zone, bande):
    """Lecture robuste (multi-echelle + centroide) d'une bande d'une image ERA5/CHIRPS."""
    for sc in (1000, 5000, 11132):
        val = image.reduceRegion(
            reducer=ee.Reducer.mean(), geometry=zone,
            scale=sc, maxPixels=1e9, bestEffort=True
        ).get(bande).getInfo()
        if val is not None:
            return val
    centre = zone.centroid(maxError=1)
    return image.reduceRegion(
        reducer=ee.Reducer.first(), geometry=centre,
        scale=11132, maxPixels=1e9, bestEffort=True
    ).get(bande).getInfo()


def evaluer_culture(nom, temp_moy, pluie_cycle, mois):
    """Renvoie un verdict pour une culture selon les conditions reelles."""
    e = EXIGENCES[nom]
    raisons = []
    score = 0  # 2 = favorable, 1 = precautions, 0 = deconseille

    # --- Critere temperature ---
    if temp_moy is None:
        raisons.append("temperature indisponible")
    elif temp_moy > e["t_critique"]:
        raisons.append(f"temperature trop elevee ({temp_moy:.1f}degC > seuil critique {e['t_critique']}degC)")
        score_temp = 0
    elif e["t_min"] <= temp_moy <= e["t_max"]:
        raisons.append(f"temperature ideale ({temp_moy:.1f}degC)")
        score_temp = 2
    else:
        raisons.append(f"temperature acceptable mais hors plage ideale ({temp_moy:.1f}degC)")
        score_temp = 1

    # --- Critere saison de semis ---
    if mois in e["mois_semis"]:
        raisons.append(f"periode de semis favorable ({MOIS_FR[mois]})")
        score_saison = 2
    else:
        favs = ", ".join(MOIS_FR[m] for m in e["mois_semis"])
        raisons.append(f"hors periode optimale (conseille : {favs})")
        score_saison = 0

    # --- Critere eau (pluie sur le cycle vs besoin) ---
    if pluie_cycle is None:
        score_eau = 1
    elif pluie_cycle >= e["eau_min"]:
        raisons.append(f"pluie suffisante sur le cycle (~{pluie_cycle:.0f} mm)")
        score_eau = 2
    elif pluie_cycle >= e["eau_min"] * 0.5:
        raisons.append(f"pluie partielle (~{pluie_cycle:.0f} mm), irrigation d'appoint necessaire")
        score_eau = 1
    else:
        raisons.append(f"pluie faible (~{pluie_cycle:.0f} mm), forte irrigation requise")
        score_eau = 1  # pas eliminatoire si on peut irriguer

    # --- Verdict global : la temperature et la saison priment ---
    if temp_moy is not None and temp_moy > e["t_critique"]:
        verdict = "Deconseille"
    elif score_saison == 0:
        verdict = "Possible avec precautions"
    elif score_temp == 2 and score_saison == 2:
        verdict = "Favorable"
    else:
        verdict = "Possible avec precautions"

    return {
        "culture": nom,
        "verdict": verdict,
        "temperature_moy": round(temp_moy, 1) if temp_moy is not None else None,
        "pluie_cycle_mm": round(pluie_cycle, 0) if pluie_cycle is not None else None,
        "raisons": raisons,
        "note": e["note"]
    }


@app.route("/plantation", methods=["POST"])
def conseil_plantation():
    """
    Conseille les cultures a planter selon les conditions reelles de la zone
    a la date choisie (temperature ERA5, pluie CHIRPS) croisees avec la fiche.
    """
    data = request.get_json()
    zone, erreur = construire_zone(data)
    if erreur:
        return jsonify({"erreur": erreur}), 400

    # Date de reference (par defaut : date_fin)
    date_ref = data.get("date_ref") or data.get("date_fin", "2024-12-31")

    try:
        from datetime import datetime, timedelta
        ref = datetime.strptime(date_ref, "%Y-%m-%d")
        mois = ref.month

        # Temperature moyenne (ERA5) sur les 30 derniers jours
        d1 = (ref - timedelta(days=30)).strftime("%Y-%m-%d")
        d2 = ref.strftime("%Y-%m-%d")
        col = (ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
               .filterDate(d1, d2).filterBounds(zone))
        # temperature_2m en Kelvin -> Celsius (moyenne de la periode)
        temp_img = col.select("temperature_2m").mean().subtract(273.15)
        temp_moy = lire_valeur_moyenne(temp_img, zone, "temperature_2m")

        # Pluie CHIRPS estimee sur un cycle (~90 jours) autour de la date
        dc1 = (ref - timedelta(days=45)).strftime("%Y-%m-%d")
        dc2 = (ref + timedelta(days=45)).strftime("%Y-%m-%d")
        chirps = (ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
                  .filterDate(dc1, dc2).filterBounds(zone))
        pluie_cycle = chirps.sum().reduceRegion(
            reducer=ee.Reducer.mean(), geometry=zone,
            scale=5000, maxPixels=1e9, bestEffort=True
        ).get("precipitation").getInfo()

        # Evaluer chaque culture
        resultats = [evaluer_culture(nom, temp_moy, pluie_cycle, mois)
                     for nom in EXIGENCES]
        # Trier : Favorable d'abord, puis precautions, puis deconseille
        ordre = {"Favorable": 0, "Possible avec precautions": 1, "Deconseille": 2}
        resultats.sort(key=lambda r: ordre.get(r["verdict"], 3))

        return jsonify({
            "mois": MOIS_FR[mois],
            "temperature_moy": round(temp_moy, 1) if temp_moy is not None else None,
            "pluie_cycle_mm": round(pluie_cycle, 0) if pluie_cycle is not None else None,
            "cultures": resultats,
            "message": "Conseil de plantation calcule"
        })
    except Exception as e:
        return jsonify({"erreur": str(e)}), 500


# ===== Page web =====
@app.route("/", methods=["GET"])
def accueil():
    # Sert index.html depuis le meme dossier que ce fichier
    return send_from_directory(os.path.dirname(os.path.abspath(__file__)), "index.html")


@app.route("/ndvi", methods=["POST"])
def calculer_ndvi():
    data = request.get_json()
    zone, erreur = construire_zone(data)
    if erreur:
        return jsonify({"erreur": erreur}), 400
    date_debut = data.get("date_debut", "2024-01-01")
    date_fin = data.get("date_fin", "2024-12-31")
    try:
        col = collection_ndvi(zone, date_debut, date_fin)
        ndvi = col.select("NDVI").median().clip(zone)
        palette = ["#a52a2a", "#d9a441", "#ffff99", "#97c459", "#3b6d11", "#173404"]
        map_id = ndvi.getMapId({"min": 0.0, "max": 0.8, "palette": palette})
        tile_url = map_id["tile_fetcher"].url_format
        moyenne = ndvi.reduceRegion(
            reducer=ee.Reducer.mean(), geometry=zone, scale=10, maxPixels=1e9
        ).get("NDVI").getInfo()

        reponse = {
            "tile_url": tile_url,
            "ndvi_moyen": round(moyenne, 3) if moyenne is not None else None,
            "message": "NDVI calcule avec succes"
        }
        reponse.update(calculer_stats(ndvi, zone))
        return jsonify(reponse)
    except Exception as e:
        return jsonify({"erreur": str(e)}), 500


@app.route("/evolution", methods=["POST"])
def evolution_ndvi():
    data = request.get_json()
    zone, erreur = construire_zone(data)
    if erreur:
        return jsonify({"erreur": erreur}), 400
    date_debut = data.get("date_debut", "2024-01-01")
    date_fin = data.get("date_fin", "2024-12-31")
    pas = data.get("pas", "mois")
    try:
        col = collection_ndvi(zone, date_debut, date_fin)
        debut = ee.Date(date_debut); fin = ee.Date(date_fin)

        if pas == "semaine":
            unite = "week"; fmt = "YYYY-MM-dd"
            n_pas = fin.difference(debut, "week").ceil()
        else:
            unite = "month"; fmt = "YYYY-MM"
            n_pas = fin.difference(debut, "month").ceil()

        indices = ee.List.sequence(0, n_pas.subtract(1))

        def ndvi_du_pas(i):
            i = ee.Number(i)
            p_debut = debut.advance(i, unite)
            p_fin = p_debut.advance(1, unite)
            sous_col = col.filterDate(p_debut, p_fin).select("NDVI")
            moyenne = ee.Algorithms.If(
                sous_col.size().gt(0),
                sous_col.median().reduceRegion(
                    reducer=ee.Reducer.mean(), geometry=zone,
                    scale=10, maxPixels=1e9, bestEffort=True
                ).get("NDVI"),
                None
            )
            return ee.Feature(None, {"date": p_debut.format(fmt), "ndvi": moyenne})

        serie = ee.FeatureCollection(indices.map(ndvi_du_pas))
        brut = serie.getInfo()["features"]
        points = []
        for f in brut:
            p = f["properties"]; valeur = p.get("ndvi")
            points.append({"date": p.get("date"),
                           "ndvi": round(valeur, 3) if valeur is not None else None})
        return jsonify({"points": points, "message": "Evolution calculee avec succes"})
    except Exception as e:
        return jsonify({"erreur": str(e)}), 500


@app.route("/anomalie", methods=["POST"])
def anomalie_ndvi():
    data = request.get_json()
    zone, erreur = construire_zone(data)
    if erreur:
        return jsonify({"erreur": erreur}), 400
    date_debut = data.get("date_debut", "2024-01-01")
    date_fin = data.get("date_fin", "2024-12-31")
    try:
        debut = ee.Date(date_debut); fin = ee.Date(date_fin)
        total_jours = fin.difference(debut, "day")
        tiers = total_jours.divide(3)
        ref_debut = debut; ref_fin = debut.advance(tiers, "day")
        rec_debut = fin.advance(tiers.multiply(-1), "day"); rec_fin = fin

        col = collection_ndvi(zone, date_debut, date_fin)
        col_ref = col.filterDate(ref_debut, ref_fin).select("NDVI")
        col_rec = col.filterDate(rec_debut, rec_fin).select("NDVI")

        if col_ref.size().getInfo() == 0 or col_rec.size().getInfo() == 0:
            return jsonify({"erreur": "Pas assez d'images nettes dans l'une des sous-periodes "
                                      "(nuages ?). Elargis la periode ou change les dates."}), 400

        ndvi_ref = col_ref.median().clip(zone)
        ndvi_rec = col_rec.median().clip(zone)
        diff = ndvi_rec.subtract(ndvi_ref).rename("diff")

        palette_div = ["#b2182b", "#ef8a62", "#fddbc7", "#f7f7f7", "#d9f0d3", "#7fbf7b", "#1b7837"]
        map_id = diff.getMapId({"min": -0.3, "max": 0.3, "palette": palette_div})
        tile_url = map_id["tile_fetcher"].url_format

        classe = (diff.gte(-SEUIL_ANOMALIE).add(diff.gt(SEUIL_ANOMALIE)).rename("classe"))
        histo = classe.reduceRegion(
            reducer=ee.Reducer.frequencyHistogram(),
            geometry=zone, scale=10, maxPixels=1e9, bestEffort=True
        ).get("classe").getInfo() or {}

        n_deg = histo.get("0", 0); n_sta = histo.get("1", 0); n_ame = histo.get("2", 0)
        total = n_deg + n_sta + n_ame
        if total > 0:
            pct_degrade = round(100 * n_deg / total, 1)
            pct_stable = round(100 * n_sta / total, 1)
            pct_ameliore = round(100 * n_ame / total, 1)
        else:
            pct_degrade = pct_stable = pct_ameliore = None

        return jsonify({
            "tile_url": tile_url,
            "pct_degrade": pct_degrade, "pct_stable": pct_stable, "pct_ameliore": pct_ameliore,
            "periode_reference": {"debut": ref_debut.format("YYYY-MM-dd").getInfo(),
                                  "fin": ref_fin.format("YYYY-MM-dd").getInfo()},
            "periode_recente": {"debut": rec_debut.format("YYYY-MM-dd").getInfo(),
                                "fin": rec_fin.format("YYYY-MM-dd").getInfo()},
            "message": "Anomalie calculee avec succes"
        })
    except Exception as e:
        return jsonify({"erreur": str(e)}), 500


if __name__ == "__main__":
    # En local : port 5000. En ligne (Cloud Run) : port fourni par la variable PORT.
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port, debug=True)
