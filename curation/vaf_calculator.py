# =============================================================================
# CALCULATEUR DE FRÉQUENCE D'ALLÈLE VARIANT (VAF)
# =============================================================================
# Ce module extrait les variants depuis les fichiers VCF et calcule
# leur Variant Allele Frequency (VAF) pour l'analyse génomique.
#
# Fonctionnalités principales:
# - Lecture des fichiers VCF (VarScan2 format)
# - Extraction des métadonnées de variants (DP, AD, FREQ)
# - Calcul du VAF (Variant Allele Frequency)
# - Export en CSV pour l'analyse ultérieure
#
# Formats supportés:
# - Fichiers VCF standards
# - Fichiers VCF VarScan2 avec métadonnées étendues

import pysam
import pandas as pd

def run():
    # Accéder directement à l'objet snakemake
    vcf_path = snakemake.input.vcf
    print(f"DEBUG: Ouverture du fichier VCF: {vcf_path}")
    
    try:
        vcf_in = pysam.VariantFile(vcf_path)
        print(f"DEBUG: Fichier VCF ouvert avec succés")
    except Exception as e:
        print(f"DEBUG: Échec de l'ouverture du VCF: {e}")
        raise e
        
    results = []
    variant_count = 0
    
    for record in vcf_in:
        variant_count += 1
        print(f"DEBUG: Traitement du variant {variant_count}: {record.chrom}:{record.pos}")
        
        # Les VCF VarScan2 contiennent généralement un échantillon
        if not record.samples:
            print(f"DEBUG: Variant {variant_count} ignoré - aucun échantillon")
            continue

        sample = record.samples[0]

        try:
            # 1. Profondeur Totale (DP)
            dp = sample.get('DP', 0)  # Total read depth at position

            # 2. Fréquence d'Allèle Variant (FREQ)
            # VarScan sort "25.4%" ou "0.254"
            freq_str = sample.get('FREQ', "0%")
            if isinstance(freq_str, str):
                vaf = float(freq_str.strip('%')) / 100.0 if '%' in freq_str else float(freq_str)
            else:
                vaf = float(freq_str)

            # 3. Profondeur Allélique (AD) - standard pour VarScan
            ad = sample.get('AD', 0)  # Read count for variant allele

            # 4. Métadonnées VarScan2 tumor-only depuis INFO
            hom = record.info.get('HOM', 0)  # Homozygous variant count
            het = record.info.get('HET', 0)  # Heterozygous variant count

            results.append({
                "CHROM": record.chrom,
                "POS": record.pos,
                "REF": record.ref,
                "ALT": record.alts[0] if record.alts else ".",
                "DP": dp,
                "AD": ad,
                "VAF": vaf,
                "HOM": hom,              # 1
                "HET": het               # 0
            })
        except Exception:
            continue

    # Créer le dataframe et exporter en CSV
    df = pd.DataFrame(results)
    print(f"DEBUG: Total des variants traités: {len(results)}")
    print(f"DEBUG: Colonnes CSV: {list(df.columns) if not df.empty else 'Aucune donnée'}")
    df.to_csv(snakemake.output.csv, index=False)
    print(f"DEBUG: CSV sauvegardé dans {snakemake.output.csv}")

if __name__ == "__main__":
    run()
