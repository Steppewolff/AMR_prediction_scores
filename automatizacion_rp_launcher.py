import automatización_rp_v2
import json
import csv

def load_scores(json_file):
    """Carga el diccionario de condiciones desde el archivo JSON."""
    with open(json_file, "r") as f:
        scores = json.load(f)
    return scores


def load_csv(csv_file):
    """Carga los registros del archivo CSV en una lista de filas."""
    records = []
    with open(csv_file, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            # Se espera que el archivo tenga 14 columnas
            records.append(row)
    return records


if __name__ == "__main__":
    # Cargar condiciones desde el JSON
    scores_json = load_scores("FuentesInformacion/SCORES_100WT.json")
    # Cargar registros del CSV (se asume que el archivo tiene 14 columnas separadas por comas)
    records = load_csv("FuentesInformacion/PA001.snps.withoutcommon.curated")

    automatización_rp_v2.main(scores_json, records)