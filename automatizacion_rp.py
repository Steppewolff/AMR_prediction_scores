import csv
import json


def effect_multiplier(effect):
    """Convierte el signo de efecto a un multiplicador (+ -> 1, - -> -1)."""
    return 1 if effect == "+" else -1


def evaluate_LOF(record):
    """
    Para condiciones LOF:
    Se comprueba si en la columna 5 (índice 4) aparece la cadena 'stop'
    (en minúsculas para comparación).
    """
    if "stop" in record[4].lower() or "indel" in record[4].lower():
        return 'LOF'
####
## Faltan por añadir los condicionales para INDELS, estos son solo para SNPs ******************************************************************************************
####

    elif 'missense_variant' in record[4].lower():
        return 'SNP'
    else:
        return 'NLOF'


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
    if "missense_variant" not in record[4].lower():
        return None
    try:
        allele_num = int(record[12])
    except ValueError:
        return None

    # Se tiene en cuenta el caso en el que no se especifican aminoácidos concretos que deban ser producto de la mutación ni la mutación es un indel o un codón STOP
    if not bool(alleles) and "stop" not in record[4].lower() and "indel" not in record[4].lower():
        return "double"  # Coincidencia: vale cualquier mutación

    # Se recorren las claves (como cadena) y se comparan como enteros, si alleles está vacío el programa no entra en el bucle
    for key, aa in alleles.items():
        try:
            ####
            ## Si lo que se registra en el Json para enumerar alelos es un rango de aminoácidos, p.ej. "80-93" : "Xxx", simplemente se comprueba si el aminoácido mutado cae dentro de ese rango, el valor Xxx no se utiliza pero se mantiene para mantener el formato del Json
            ####
            if '-' in key:
                keys = key.split('-')
                if int(keys[0]) <= allele_num <= int(keys[1]):
                    return "double"  # Coincidencia: se usará el valor doble
            elif int(key) == allele_num:
                values = aa if isinstance(aa, list) else [aa]
                for value in values:
                ####
                ## Hay que revisar si los casos del siguiente condicional leen las columnas que tienen que leer y son suficientes para
                ## identificar correctamente los casos que aparecen en las condiciones de Scores ******************************************************************************************
                ####
                    # Se compara el alelo esperado con las 3 últimas letras de la anotación proteica
                    if value == 'Xxx' and "stop" not in record[4].lower() and "indel" not in record[4].lower():
                        return "double"  # Coincidencia: vale cualquier mutación
                    elif record[9][-3:] == value:
                        return "double"  # Coincidencia: se usará el valor doble

                return "simple"
        except ValueError:
            continue

    ####
    ## ¿Duplicado con el anterior return "simple" si pasa por uno podría pasar por el otro?
    ####
    # return "simple"


def evaluate_condition(record, condition):
    """
    Evalúa la condición (sin considerar reguladores) para un registro dado.
    Devuelve el score (con signo) a sumar para ese registro según la condición.
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

    elif mut_type == "LOFN":
        if evaluate_LOF == 'LOF':
            reg_result['reg_eval'] = 0
            reg_result['value'] = 0
        ## Si es un SNP en los genes de bomba de mex_s, se puede considerar directamente que está activa? o hay que evaluar algo más?  ******************************************************************************************
        elif evaluate_LOF == 'SNP':
            reg_result['reg_eval'] = 1
            reg_result['value'] = eff * condition["simple_value"]
        else:
            reg_result['reg_eval'] = 1
            reg_result['value'] = eff * condition["doble_value"]

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
        if regulator_pass:
            break
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
        # Se suma tanto el score del nivel principal como un "bonus" (usando el valor doble)
        # return main_score + effect_multiplier(condition["effect"]) * condition["doble_value"]
        return main_score + regulator_value
    elif not regulator_pass:
        # Si no se cumple el nivel regulador, se suma el score principal
        # return main_score if main_score != 0 else effect_multiplier(condition["effect"]) * condition["simple_value"]
        return main_score

    ###
    # Este caso es necesario?*************************************************************************************************
    ###
    elif regulator_pass and main_score == 0:
        return effect_multiplier(condition["effect"]) * condition["simple_value"]
    # Caso por defecto
    return main_score


def main(scores_json, records):
    # # Cargar condiciones desde el JSON
    # scores_json = load_scores("FuentesInformacion/SCORES_100WT.json")
    # # Cargar registros del CSV (se asume que el archivo tiene 14 columnas separadas por comas)
    # records = load_csv("FuentesInformacion/PA001.snps.withoutcommon.curated")

    # Inicializar diccionario de scores para cada antibiótico
    antibiotics = ["CIP", "CAZ", "MER", "C/T", "TOB"]
    final_scores = {ab: 0 for ab in antibiotics}
    final_score_eval = {ab: [] for ab in antibiotics}

    gene_records = [r for r in records if r[7] in scores_json]
    gene_list = [r[7] for r in gene_records]

    # Recorrer cada gen (locus) definido en el JSON
    for gene, conditions in scores_json.items():
        # Caso en el que son varios genes principales, no los reguladores los que determinan si se suma score o no
        if ',' in gene:
            genes = gene.split(',')
            gene_records_test = [r for r in records if r[7] in genes]
            gene_active = True
            score = 0

            # Siempre van a tener reguladores
            for gene_record in gene_records_test:
                if evaluate_LOF(gene_record) == 'LOF':
                    gene_active = False
                    break

            if gene_active:
                for ab, cond in conditions.items():
                    for sub_gene in genes:
                        score = evaluate_regulators(sub_gene, cond, records)
                        final_scores[ab] += score
                        final_score_eval[ab].append({sub_gene: score})

        # Para cada antibiótico (condición) en el gen
        if gene in gene_list:
            mutations = [r for r in gene_records if r[7] == gene]
            for ab, cond in conditions.items():
                # Se verifica el tipo de evaluación según el campo "regulators"
                if cond.get("regulators") == "NO":
                    # Se filtran los registros cuya columna 8 (índice 7) coincide con el gen
                    for mutation in mutations:
                        score = evaluate_condition(mutation, cond)['value']
                        final_scores[ab] += score
                        final_score_eval[ab].append({gene: score})

                elif cond.get("regulators") == "YES":
                    # Para condiciones con reguladores, se evalúa con la función especializada
                    score = evaluate_regulators(gene, cond, records)
                    final_scores[ab] += score
                    final_score_eval[ab].append({gene: score})


    # Mostrar el score final por antibiótico
    print("Total scores by ATB:")
    for ab, score in final_scores.items():
        print(f"{ab}: {score}")
    for ab, scores in final_score_eval.items():
        print(f"{ab}: {scores}")