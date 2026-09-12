# passerelle.py — Passerelle AgriScan
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
