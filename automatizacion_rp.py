import pandas as pd
import re


def parse_single_value(value_str):
    """
    Convierte un valor numérico (en formato string) eliminando posibles signos y
    sustituyendo la coma decimal por punto.
    """
    value_clean = value_str.strip()
    if value_clean and value_clean[0] in ['+', '-']:
        value_clean = value_clean[1:]
    value_clean = value_clean.replace(',', '.')
    try:
        return float(value_clean)
    except Exception as e:
        raise ValueError(f"No se pudo convertir el valor '{value_str}': {e}")


def parse_doble_score_value(value_str, observations, snp_aa):
    """
    Cuando VALUE contiene dos valores separados por '/' y se ha indicado DOBLE SCORE,
    se selecciona el valor a usar:
      - Se extraen los dos valores.
      - Si OBSERVACIONES incluye una lista de mutaciones (p.ej. DOBLE SCORE (L14Q, V101M, ...)),
        se extraen los números de aminoácido y si el número de AA_POS (primer número de la cadena)
        se encuentra en esa lista se usa el primer valor (mayor); si no, se usa el segundo (menor).
      - Si no se encuentra listado, se usa por defecto el primer valor.
    """
    parts = value_str.split('/')
    if len(parts) != 2:
        return parse_single_value(value_str)
    try:
        v1 = parse_single_value(parts[0])
        v2 = parse_single_value(parts[1])
    except Exception as e:
        raise ValueError(f"Error al parsear doble score en '{value_str}': {e}")

    # Buscar listado de mutaciones tras "DOBLE SCORE"
    match = re.search(r'DOBLE SCORE\s*\((.*?)\)', observations, re.IGNORECASE)
    if match:
        mutations_str = match.group(1)  # Ejemplo: "L14Q, V101M, L137P, A138T, A168V, Q232E, G361R"
        mutation_positions = []
        for num in re.findall(r'\d+', mutations_str):
            try:
                mutation_positions.append(int(num))
            except:
                continue
        if mutation_positions and snp_aa in mutation_positions:
            return v1  # Se usa el primer valor (mayor)
        else:
            return v2  # Se usa el segundo valor (menor)
    else:
        # Si no se incluye listado, se usa siempre el primer valor.
        return v1


def check_region_condition(observations, snp_aa):
    """
    Verifica si en OBSERVACIONES aparece SOLO REGION QRDR o SOLO QRDR y
    se especifica un intervalo (por ejemplo, "466-468") o una lista de posiciones
    (por ejemplo, "(81, 83, 86, 87, 106)"). Devuelve True si el valor de snp_aa cumple la condición.
    """
    # Buscar patrón de intervalo: "466-468"
    range_match = re.search(r'(\d+)\s*-\s*(\d+)', observations)
    if range_match:
        low = int(range_match.group(1))
        high = int(range_match.group(2))
        return low <= snp_aa <= high

    # Buscar patrón de lista: "(81, 83, 86, 87, 106)"
    list_match = re.search(r'\(([\d,\s]+)\)', observations)
    if list_match:
        nums_str = list_match.group(1)
        try:
            numbers = [int(x.strip()) for x in nums_str.split(',') if x.strip().isdigit()]
        except:
            numbers = []
        return snp_aa in numbers

    return False


def check_stop_condition(observations, snp_effect):
    """
    Verifica si en OBSERVACIONES aparece SOLO STOP CODON/FS y si en la columna EFFECT del SNP
    se encuentra la cadena 'missense_variant'.
    """
    if "SOLO STOP CODON/FS" in observations.upper():
        return "missense_variant" in snp_effect.lower()
    return False


def evaluate_scores(scores_file, snps_file):
    # Cargar la pestaña "SCORES" del archivo de Excel
    df_scores = pd.read_excel(scores_file, sheet_name="SCORES")
    # Cargar el archivo CSV con los SNPs
    df_snps = pd.read_csv(snps_file)

    # Inicializar diccionario para acumular scores por antibiótico
    antibiotics = ["CAZ", "MER", "CIP", "C/T", "TOB"]
    antibiotic_scores = {ab: 0 for ab in antibiotics}

    # Recorrer cada condición (registro) en el archivo de scores
    for idx, score_row in df_scores.iterrows():
        locus = score_row["LOCUS"]
        antibiotic = score_row["ANTIBIOTIC"]
        value_str = str(score_row["VALUE"]).strip()
        effect_str = str(score_row["EFFECT"]).strip()
        observations = str(score_row["OBSERVACIONES"]).strip() if not pd.isna(score_row["OBSERVACIONES"]) else ""

        # Determinar el multiplicador a partir de EFFECT (se ignora el signo en VALUE)
        multiplier = 1
        if '-' in effect_str:
            multiplier = -1
        elif '+' in effect_str:
            multiplier = 1

        # Buscar en snps.csv aquellos registros cuyo LOCUS_TAG coincide con el LOCUS
        matching_snps = df_snps[df_snps["LOCUS_TAG"] == locus]

        # Variable para controlar si la condición se ha cumplido para este score
        condition_met = False
        chosen_value = None

        # Iterar sobre los registros del SNP para el locus en cuestión
        for i, snp_row in matching_snps.iterrows():
            # Extraer la posición de aminoácido (se espera formato "num_total", por ejemplo "123/456")
            aa_pos_field = str(snp_row["AA_POS"]).strip()
            aa_parts = aa_pos_field.split('/')
            try:
                snp_aa = int(aa_parts[0].strip())
            except:
                continue  # Si no se puede parsear, se omite el registro

            # Evaluar condiciones adicionales según OBSERVACIONES
            if "DOBLE SCORE" in observations.upper():
                if "/" in value_str:
                    chosen_value = parse_doble_score_value(value_str, observations, snp_aa)
                else:
                    chosen_value = parse_single_value(value_str)
                condition_met = True
            elif ("SOLO REGION QRDR" in observations.upper()) or ("SOLO QRDR" in observations.upper()):
                if check_region_condition(observations, snp_aa):
                    # Si por error VALUE tiene dos valores, se toma el primero
                    if "/" in value_str:
                        chosen_value = parse_single_value(value_str.split('/')[0])
                    else:
                        chosen_value = parse_single_value(value_str)
                    condition_met = True
                else:
                    condition_met = False
            elif "SOLO STOP CODON/FS" in observations.upper():
                if check_stop_condition(observations, str(snp_row["EFFECT"])):
                    if "/" in value_str:
                        chosen_value = parse_single_value(value_str.split('/')[0])
                    else:
                        chosen_value = parse_single_value(value_str)
                    condition_met = True
                else:
                    condition_met = False
            else:
                # Sin condiciones adicionales: se asume que basta con que el LOCUS coincida.
                if "/" in value_str:
                    # Si aparecen dos valores y no se indicó DOBLE SCORE, se toma el primero por defecto.
                    chosen_value = parse_single_value(value_str.split('/')[0])
                else:
                    chosen_value = parse_single_value(value_str)
                condition_met = True

            # Si se cumple la condición para este registro SNP, se actualiza el score para el antibiótico
            if condition_met and chosen_value is not None:
                if antibiotic in antibiotic_scores:
                    antibiotic_scores[antibiotic] += multiplier * chosen_value
                else:
                    antibiotic_scores[antibiotic] = multiplier * chosen_value
                # Se cuenta solo una vez por cada registro de scores; se rompe el bucle de SNPs.
                break

    return antibiotic_scores


if __name__ == '__main__':
    scores_file = 'SCORES_100WT.xlsx'
    snps_file = 'snps.csv'

    results = evaluate_scores(scores_file, snps_file)
    print("Resultados por antibiótico:")
    for ab, score in results.items():
        print(f"{ab}: {score}")