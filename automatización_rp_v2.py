import csv
import json


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


def effect_multiplier(effect):
    """Convierte el signo de efecto a un multiplicador (+ -> 1, - -> -1)."""
    return 1 if effect == "+" else -1


def evaluate_LOF(record):
    """
    Para condiciones LOF:
    Se comprueba si en la columna 5 (índice 4) aparece la cadena 'stop'
    (en minúsculas para comparación).
    """
    if "stop" in record[4].lower():
        return 'LOF'
####
## Faltan por añadir los condicionales para INDELS, estos son solo para SNPs ******************************************************************************************
####

    if 'missense_Variant' in record[4].lower():
        return 'SNP'


def evaluate_GOF(record, alleles):
    """
    Para condiciones GOF/GOFO:
    - Se verifica que en la columna 5 (índice 4) aparezca 'missense_variant'.
    - Se extrae el valor numérico de la columna 13 (índice 12) y se compara con las claves
      de 'alleles' (convertidas a entero).
    - Si hay coincidencia, se compara el valor esperado (cadena de 3 caracteres) con las
      3 últimas letras de la anotación proteica (se asume en la columna 10, índice 9).
    Devuelve:
      - "double" si se cumple la condición de coincidencia de alelo.
      - "simple" si se cumple la condición de variante missense pero no coincide el alelo.
      - None si no se cumple.
    """
    if "missense_variant" not in record[4]:
        return None
    try:
        allele_num = int(record[12])
    except ValueError:
        return None
    # Se recorren las claves (como cadena) y se comparan como enteros
    for key, aa in alleles.items():
        try:
            if int(key) == allele_num:
                # Se compara el alelo esperado con las 3 últimas letras de la anotación proteica
                if record[9][-3:] == aa:
                    return "double"  # Coincidencia: se usará el valor doble
                else:
                    return "simple"
        except ValueError:
            continue
    return "simple"


def evaluate_condition(record, condition):
    """
    Evalúa la condición (sin considerar reguladores) para un registro dado.
    Devuelve el puntaje (con signo) a sumar para ese registro según la condición.
    Si no se cumple, devuelve 0.
    """
    mut_type = condition["mutation_type"]
    eff = effect_multiplier(condition["effect"])
    reg_result = {'reg_eval' : 0, 'value' : 0}

    if mut_type == "LOF":
        if evaluate_LOF == 'LOF':
            reg_result['reg_eval'] = 1
            reg_result['value'] = eff * condition["doble_value"]
        elif evaluate_LOF == 'SNP':
            reg_result['reg_eval'] = 1
            reg_result['value'] = eff * condition["simple_value"]
        else:
            reg_result['reg_eval'] = 0
            reg_result['value'] = 0

    elif mut_type in ("GOF", "GOFO"):
        result = evaluate_GOF(record, condition["alleles"])
        if result is None:
            reg_result['reg_eval'] = 0
            reg_result['value'] = 0
        if result == "double":
            reg_result['reg_eval'] = 1
            reg_result['value'] = eff * condition["doble_value"]
        else:
            # Para GOF se suma simple_value si no coincide el alelo;
            # Para GOFO, si no hay coincidencia, no se suma nada.
            if mut_type == "GOF":
                reg_result['reg_eval'] = 1
                reg_result['value'] = eff * condition["simple_value"]
            else:
                reg_result['reg_eval'] = 0
                reg_result['value'] = 0
    else:
        print(f'Error in mut_type value for {record[7]} it is none of LOF, GOF or GOFO')

    return reg_result


def evaluate_regulators(gene, condition, records):
    """
    Evalúa la condición cuando 'regulators' es "YES".
    Se consideran dos niveles:
      1. Nivel de los reguladores (sub-loci): se revisa si existe al menos un registro
         cuyo locus (columna 8) coincide con alguna clave en 'loci' y que cumpla la condición
         del sublocus.
      2. Nivel superior: se evalúa la condición principal en los registros donde el locus (columna 8)
         coincide con el gen.
    La combinación de ambos niveles se implementa (de forma simplificada) de la siguiente manera:
      - Si existe un registro en algún sublocus que cumpla la condición (regulator_pass)
        y además se cumple la condición principal (main_score ≠ 0), se suma:
             main_score + (bonus de doble_value de la condición principal)
      - Si no se cumple el primer nivel, se suma el main_score (o el simple_value de la condición principal).
      - Si se cumple el nivel de reguladores pero no la condición principal, se suma el simple_value.
    (La lógica puede ajustarse según la interpretación de los requisitos.)
    """
    # Registros correspondientes al gen principal
    main_records = [r for r in records if r[7] == gene]
    # Se obtienen las claves de los sub-loci
    subloci = condition.get("loci", {})
    # Se recopilan todos los registros que pertenezcan a alguno de los sub-loci
    sub_records = []
    for sub in subloci:
        sub_records.extend([r for r in records if r[7] == sub])

    # Evaluar los sub-loci: verificar si existe al menos un registro que cumpla la condición del sublocus y sumar su valor si se cumplen las condiciones
    regulator_pass = False
    regulator_value = 0
    for sub in sub_records:
    # for sub in subloci:
        sub_cond = subloci[sub[7]]
        sub_recs = [r for r in records if r[7] == sub[7]]
        for rec in sub_recs:
####
## En la siguiente condición solamente se tiene en cuenta el primer regulador cuyo efecto es diferente de 0, si hay 2 con efecto deberían sumarse los 2?? ************************************************************
####
            regulator_condition = evaluate_condition(rec, sub_cond)
            if regulator_condition['reg_eval'] != 0:
                regulator_pass = True
                regulator_value += regulator_condition['value']

                break
        if regulator_pass:
            break

    # Evaluar la condición principal en los registros del gen
    main_score = 0

###
# En un locus concreto, si hay más de una mutación ¿como se evalúa? Puede haber una con un efecto y otra con otro, ¿se suman ambas? ¿Se escoge la más fuerte?
###
    for rec in main_records:
        main_score += evaluate_condition(rec, condition)['value']

    # Aplicar una lógica combinada (según la descripción)
    # if regulator_pass and main_score != 0:
    if regulator_pass:
        # Se suma tanto el puntaje del nivel principal como un "bonus" (usando el valor doble)
        # return main_score + effect_multiplier(condition["effect"]) * condition["doble_value"]
        return main_score + regulator_value
    elif not regulator_pass:
        # Si no se cumple el nivel regulador, se suma el puntaje principal
        # return main_score if main_score != 0 else effect_multiplier(condition["effect"]) * condition["simple_value"]
        return main_score
    elif regulator_pass and main_score == 0:
        return effect_multiplier(condition["effect"]) * condition["simple_value"]
    # Caso por defecto
    return main_score


def main(scores_json, records):
    # # Cargar condiciones desde el JSON
    # scores_json = load_scores("FuentesInformacion/SCORES_100WT.json")
    # # Cargar registros del CSV (se asume que el archivo tiene 14 columnas separadas por comas)
    # records = load_csv("FuentesInformacion/PA001.snps.withoutcommon.curated")

    # Inicializar diccionario de puntajes para cada antibiótico
    antibiotics = ["CIP", "CAZ", "MER", "C/T", "TOB"]
    final_scores = {ab: 0 for ab in antibiotics}

    # Recorrer cada gen (locus) definido en el JSON
    for gene, conditions in scores_json.items():
        # Para cada antibiótico (condición) en el gen
        for ab, cond in conditions.items():
            # Se verifica el tipo de evaluación según el campo "regulators"
            # if cond.get("regulators", "NO") == "NO":
            if cond.get("regulators") == "NO":
                # Se filtran los registros cuya columna 8 (índice 7) coincide con el gen
                gene_records = [r for r in records if r[7] == gene]
                for rec in gene_records:
                    final_scores[ab] += evaluate_condition(rec, cond)['value']
            elif cond.get("regulators") == "YES":
                # Para condiciones con reguladores, se evalúa con la función especializada
                score = evaluate_regulators(gene, cond, records)
                final_scores[ab] += score
                pass

            pass

    # Mostrar el puntaje final por antibiótico
    print("Total scores by ATB:")
    for ab, score in final_scores.items():
        print(f"{ab}: {score}")