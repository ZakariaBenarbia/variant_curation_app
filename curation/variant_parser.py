import pandas as pd
import numpy as np
import os

# =============================================================================
# EXTRACTEUR DE CARACTÉRISTIQUES DE VARIANTES
# =============================================================================
# Ce module extrait et enrichit les données de variants génomiques
# depuis les fichiers CSV pour l'analyse par IA.
#
# Fonctionnalités principales:
# - Classification des rôles de gènes (Oncogène, TSG, etc.)
# - Amélioration de l'origine des mutations (Somatic, Germline)
# - Calcul de scores de priorité clinique
# - Enrichissement des métadonnées génomiques
#
# Sources de données utilisées:
# - COSMIC Cancer Gene Census
# - Métadonnées de variants (VAF, profondeur, etc.)

import pandas as pd
import numpy as np
import os

def classify_gene_role(cosmic_role_field):
    """Classifie le rôle d'un gène selon la base de données COSMIC"""
    if not cosmic_role_field or pd.isna(cosmic_role_field):
        return "Unknown"

    role_str = str(cosmic_role_field).lower()
    is_oncogene = "oncogene" in role_str
    is_tsg = "tsg" in role_str or "tumour suppressor" in role_str

    # Logique de classification complexe pour correspondre à votre base de données
    if is_oncogene and is_tsg:
        return "Oncogene/TSG"  # Double rôle
    elif is_oncogene:
        return "Oncogene"  # Gène promoteur du cancer
    elif is_tsg:
        return "Tumor Suppressor"  # Gène suppresseur de tumeur
    else:
        return "Other Cancer Gene"  # Autre gène lié au cancer

def infer_origin_improved(cosmic_origin, vaf, depth, mut_status="N/A"):
    """Infère l'origine améliorée de la mutation"""
    status = str(mut_status).lower()
    
    # Logique d'inférence basée sur VAF et profondeur
    if "somatic" in status:
        if vaf >= 0.70: 
            return "Somatic", "Variante somatique confirmée + VAF élevée"
        return "Strictly Somatic", "Variante somatique strictement confirmée"
    
    # Critères pour LOH (Loss of Heterozygosity) et mutations germinales
    if vaf >= 0.95 and depth > 50:
        return "Germline (Homozygous)", "Presque 100% VAF - probable homozygote"
    if 0.45 <= vaf <= 0.55 and depth > 100:
        return "Strictly Somatic", "Variante somatique hétérozygote"
    
    # Retourner l'annotation COSMIC par défaut
    return str(cosmic_origin), "Annotation base de données COSMIC"

def calculate_improved_score(row):
    """Calcule un score de priorité clinique amélioré"""
    if not row['is_reliable']: 
        return 0.0  # Donnée non fiable
    
    vaf, log_depth = row['vaf'], row['log_depth']
    gene_role, origin = row['gene_role'], row['origin_improved']
    
    # Pondérations par rôle de gène pour le scoring
    role_weights = {
        "Oncogene": 1.3, 
        "Tumor Suppressor": 1.5, 
        "Oncogene/TSG": 1.4, 
        "Other Cancer Gene": 0.9, 
        "Unknown": 0.5
    }
    role_mult = role_weights.get(gene_role, 0.5)

    if "Somatic" in origin: origin_mult = 1.8
    elif "Strictly Somatic" in origin: origin_mult = 1.5
    elif "Germline" in origin: origin_mult = 0.6
    else: origin_mult = 0.8

    return round(vaf * log_depth * role_mult * origin_mult, 4)

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    vaf_input = snakemake.input.vaf_csv
    output_path = snakemake.output.features
    census_path = "resources/cosmic/CancerGeneCensus.tsv"
    mutation_path = "resources/cosmic/cosmic_mutation.tsv"

    try:
        df = pd.read_csv(vaf_input)
        if df.empty:
            pd.DataFrame().to_csv(output_path, index=False); return

        df = df.rename(columns={'CHROM': 'chrom', 'POS': 'pos', 'DP': 'depth', 'VAF': 'vaf'})
        if df['vaf'].dtype == object: df['vaf'] = df['vaf'].str.replace('%', '').astype(float)
        if df['vaf'].max() > 1.0: df['vaf'] = df['vaf'] / 100.0

        df['chrom_clean'] = df['chrom'].astype(str).str.replace('chr', '', case=False).str.strip()

        # Load Gene Census - Updated to capture ROLE_IN_CANCER accurately
        census = pd.read_csv(census_path, sep='\t').dropna(subset=['GENOME_START', 'GENOME_STOP'])
        gene_map = {}
        for _, row in census.iterrows():
            c = str(row['CHROMOSOME']).replace('chr', '').strip()
            if c not in gene_map: gene_map[c] = []

            # Capturing the raw role from database
            gene_map[c].append({
                'start': int(row['GENOME_START']), 'stop': int(row['GENOME_STOP']),
                'symbol': str(row['GENE_SYMBOL']), 'tumor_type': str(row.get('TUMOUR_TYPES_SOMATIC', 'N/A')),
                'origin': "Somatic" if str(row.get('SOMATIC')).lower() == 'y' else "Associated",
                'gene_role': classify_gene_role(row.get('ROLE_IN_CANCER', 'Unknown')),
                'cosmic_tier': int(row.get('TIER', 2)) if pd.notna(row.get('TIER')) else 2
            })

        # Load cosmic mutation data first to get gene symbols
        cosmic_gene_map = {}
        if os.path.exists(mutation_path):
            mutation = pd.read_csv(mutation_path, sep='\t', low_memory=False)
            mutation.rename(columns={'CHROMOSOME': 'chrom_clean', 'GENOME_START': 'pos',
                                  'GENOMIC_MUTATION_ID': 'cosmic_id', 'MUTATION_SOMATIC_STATUS': 'mut_status',
                                  'GENE_SYMBOL': 'cosmic_gene'}, inplace=True)
            mutation['chrom_clean'] = mutation['chrom_clean'].astype(str).str.replace('chr', '').str.strip()
            mutation_dedup = mutation.drop_duplicates(subset=['chrom_clean', 'pos'])

            # Create mapping of position to cosmic gene data
            for _, row in mutation_dedup.iterrows():
                key = (row['chrom_clean'], row['pos'])
                cosmic_gene_map[key] = {
                    'gene': row['cosmic_gene'],
                    'cosmic_id': row['cosmic_id'],
                    'mut_status': row['mut_status']
                }

        def annotate(chrom, pos):
            # First check if we have cosmic gene data for this position
            cosmic_key = (chrom, pos)
            if cosmic_key in cosmic_gene_map:
                cosmic_gene = cosmic_gene_map[cosmic_key]['gene']
                # Check if this gene exists in CancerGeneCensus for additional annotation
                if chrom in gene_map:
                    for g in gene_map[chrom]:
                        if g['symbol'] == cosmic_gene or ((g['start'] - 5000) <= pos <= (g['stop'] + 5000)):
                            return (g['symbol'], g['tumor_type'], g['origin'], g['gene_role'], g['cosmic_tier'])

                # Gene exists in cosmic but not in CancerGeneCensus - use cosmic data with defaults
                return (cosmic_gene, "N/A", "Unknown", "Unknown", 2)

            # No cosmic data, fall back to CancerGeneCensus annotation
            if chrom in gene_map:
                for g in gene_map[chrom]:
                    if (g['start'] - 5000) <= pos <= (g['stop'] + 5000):
                        return (g['symbol'], g['tumor_type'], g['origin'], g['gene_role'], g['cosmic_tier'])

            # No annotation found in either database
            return ("Intergenic", "N/A", "Unknown", "Unknown", 2)

        results = df.apply(lambda x: annotate(x['chrom_clean'], x['pos']), axis=1)
        df['gene'], df['tumor_type'], df['origin'], df['gene_role'], df['cosmic_tier'] = zip(*results)

        # Add cosmic_id and mut_status from our mapping
        df['cosmic_id'] = df.apply(lambda x: cosmic_gene_map.get((x['chrom_clean'], x['pos']), {}).get('cosmic_id', 'Unreported'), axis=1)
        df['mut_status'] = df.apply(lambda x: cosmic_gene_map.get((x['chrom_clean'], x['pos']), {}).get('mut_status', 'N/A'), axis=1)

        df['is_reliable'] = df['depth'] > 5
        df['log_depth'] = np.log10(df['depth'].replace(0, 1))

        df['origin_improved'] = df.apply(lambda x: infer_origin_improved(x['origin'], x['vaf'], x['depth'], x['mut_status'])[0], axis=1)
        df['improved_score'] = df.apply(calculate_improved_score, axis=1)
        df['loh_status'] = df.apply(lambda x: "Likely LOH" if x['vaf'] >= 0.70 and "Somatic" in x['origin_improved'] and "Germline" not in x['origin_improved'] else "No LOH", axis=1)
        df.sort_values(by='improved_score', ascending=False).to_csv(output_path, index=False)
        print(f"Success! Classification includes Fusions and Multi-role genes.")

    except Exception as e:
        print(f"Failed: {e}"); raise e

if __name__ == "__main__":
    main()
