# ============================
# INFORME MASIVO DE GENOTIPADO KIR DESDE KALLISTO
# ============================

import pandas as pd   # Librería principal para manejo de tablas
import os             # Para gestionar rutas y recorrer directorios

# ------------------------------------------------
# CONFIGURACIÓN DE RUTAS Y PARÁMETROS
# ------------------------------------------------

# Ruta base del proyecto
base_path = "/home/davidlazaro/Desktop/hackaton_lamarato2025"

# Carpeta donde kallisto guarda una subcarpeta por muestra
results_dir = os.path.join(base_path, "results/kallisto")

# Fichero que traduce target_id -> allele_name (diccionario KIR)
desc_file = os.path.join(base_path, "db/kir/KIR_descriptions.txt")

# Umbral mínimo de evidencia para descartar ruido
umbral_ruido = 100

# ------------------------------------------------
# 1) CARGA DEL DICCIONARIO DE ALELOS (UNA SOLA VEZ)
# ------------------------------------------------

# Leemos el fichero de descripciones de alelos KIR
desc = pd.read_csv(
    desc_file,
    sep=r'\s+',                   # Separado por espacios
    header=None,                  # No tiene cabecera
    names=['target_id', 'allele_name'],
    usecols=[0, 1],               # Solo target_id y allele_name
    engine='python'
)

# Normalizamos target_id para evitar errores en el merge
desc['target_id'] = desc['target_id'].astype(str).str.strip()

# ------------------------------------------------
# 2) FUNCIÓN PRINCIPAL: PROCESAR UNA MUESTRA
# ------------------------------------------------

def procesar_muestra(nombre_muestra, ruta_tsv):
    """
    Procesa un archivo abundance.tsv de kallisto y devuelve
    una tabla resumen con:
      - Gen KIR
      - Tipo de cigocidad (HOMO / HETERO)
      - Genotipo inferido
    """

    # Leemos la cuantificación de kallisto
    df = pd.read_csv(ruta_tsv, sep='\t')

    # Normalizamos target_id
    df['target_id'] = df['target_id'].astype(str).str.strip()
    
    # Cruzamos con el diccionario para obtener el nombre del alelo
    result = pd.merge(df, desc, on='target_id')

    # Extraemos el gen a partir del nombre del alelo
    # Ejemplo: KIR3DL1*0010101 -> KIR3DL1
    result['gene'] = result['allele_name'].apply(
        lambda x: x.split('*')[0]
    )
    
    # ------------------------------------------------
    # FILTRADO DE RUIDO
    # ------------------------------------------------

    # Eliminamos alelos con poca evidencia (est_counts bajos)
    df_clean = result[result['est_counts'] > umbral_ruido].copy()
    
    # Lista donde acumularemos el resumen por gen
    resumen_paciente = []

    # ------------------------------------------------
    # INFERENCIA DE GENOTIPO POR GEN KIR
    # ------------------------------------------------

    for gene, group in df_clean.groupby('gene'):

        # Ordenamos alelos por número de lecturas estimadas
        group = group.sort_values(by='est_counts', ascending=False)

        # Alelo principal (mayor evidencia)
        a1 = group.iloc[0]['allele_name']
        c1 = group.iloc[0]['est_counts']
        
        # Regla de cigocidad:
        # Si el segundo alelo tiene >25% de la señal del primero,
        # se considera heterocigoto
        if len(group) > 1 and (group.iloc[1]['est_counts'] / c1) > 0.25:
            genotipo = f"{a1} / {group.iloc[1]['allele_name']}"
            cigocidad = "HETERO"
        else:
            genotipo = a1
            cigocidad = "HOMO"
        
        # Añadimos el resultado para este gen
        resumen_paciente.append([gene, cigocidad, genotipo])
    
    # Devolvemos una tabla limpia y legible
    return pd.DataFrame(
        resumen_paciente,
        columns=['GEN', 'TIPO', 'GENOTIPO']
    )

# ------------------------------------------------
# 3) BUCLE PRINCIPAL: INFORME PARA TODAS LAS MUESTRAS
# ------------------------------------------------

print("\n" + "="*60)
print("  SISTEMA DE ANÁLISIS GENÉTICO KIR - INFORME MASIVO")
print("="*60)

# Listamos todas las carpetas de muestras
muestras_encontradas = sorted(os.listdir(results_dir))

for carpeta in muestras_encontradas:

    # Ruta al abundance.tsv de la muestra
    ruta_tsv = os.path.join(results_dir, carpeta, "abundance.tsv")
    
    # Solo procesamos si existe el archivo
    if os.path.exists(ruta_tsv):

        # Cabecera visual por muestra
        print(f"\n»»» MUESTRA: {carpeta}")
        print("-" * 50)
        
        # Procesamos la muestra
        info_genotipo = procesar_muestra(carpeta, ruta_tsv)
        
        if not info_genotipo.empty:
            # Imprimimos la tabla sin el índice de pandas
            print(info_genotipo.to_string(index=False))
        else:
            # Caso sin genes por encima del umbral
            print("AVISO: No se detectaron genes por encima del umbral de ruido.")
            
        print("-" * 50)

# ------------------------------------------------
# FIN DEL PROCESAMIENTO
# ------------------------------------------------

print("\n" + "="*60)
print("  FIN DEL PROCESAMIENTO")
print("="*60)
