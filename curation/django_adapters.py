"""
Adapter functions to integrate the heme-variant-pipeline scripts with Django
"""
import pandas as pd
import os
from .vaf_calculator import run as calculate_vaf
from .variant_parser import main as extract_features
from .ml_classifier import main as score_variants
from .report_generator import main as generate_clinical_report


class MockSnakemake:
    """Mock snakemake object for script compatibility"""
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)


def parse_vcf_to_csv(vcf_path, csv_path):
    """Convert VCF to CSV using the vaf_calculator script"""
    # Create mock snakemake object
    mock_snakemake = MockSnakemake(
        input=MockSnakemake(vcf=vcf_path),
        output=MockSnakemake(csv=csv_path)
    )
    
    # Temporarily replace snakemake import
    import sys
    import builtins
    original_snakemake = getattr(sys.modules.get('vaf_calculator', {}), 'snakemake', None)
    
    # Add snakemake to the module's namespace
    import importlib
    vaf_module = importlib.import_module('.vaf_calculator', package='curation')
    vaf_module.snakemake = mock_snakemake
    
    try:
        calculate_vaf()
        return True
    except Exception as e:
        print(f"Error in VCF to CSV conversion: {e}")
        return False


def process_variant_features(csv_path, features_path):
    """Extract features from CSV using the variant_parser script"""
    mock_snakemake = MockSnakemake(
        input=MockSnakemake(vaf_csv=csv_path),
        output=MockSnakemake(features=features_path)
    )
    
    import importlib
    parser_module = importlib.import_module('.variant_parser', package='curation')
    parser_module.snakemake = mock_snakemake
    
    try:
        extract_features()
        return True
    except Exception as e:
        print(f"Error in feature extraction: {e}")
        return False


def apply_ai_scoring(features_path, scores_path):
    """Apply AI scoring using the ml_classifier script"""
    mock_snakemake = MockSnakemake(
        input=MockSnakemake(features=features_path),
        output=MockSnakemake(scores=scores_path)
    )
    
    import importlib
    ml_module = importlib.import_module('.ml_classifier', package='curation')
    ml_module.snakemake = mock_snakemake
    
    try:
        score_variants()
        return True
    except Exception as e:
        print(f"Error in AI scoring: {e}")
        return False


def generate_report(scores_path, report_path, sample_id="sample"):
    """Generate HTML report using the report_generator script"""
    # Create output paths for charts
    vaf_plot = report_path.replace('.html', '_vaf_dist.png')
    gene_role_chart = report_path.replace('.html', '_gene_roles.png')
    loh_pie_chart = report_path.replace('.html', '_mutation_origin.png')
    
    mock_snakemake = MockSnakemake(
        input=MockSnakemake(scores=scores_path),
        output=MockSnakemake(
            report=report_path,
            vaf_plot=vaf_plot,
            gene_role_bar_chart=gene_role_chart,
            loh_pie_chart=loh_pie_chart
        ),
        wildcards=MockSnakemake(sample=sample_id)
    )
    
    import importlib
    report_module = importlib.import_module('.report_generator', package='curation')
    report_module.snakemake = mock_snakemake
    
    try:
        generate_clinical_report()
        return True
    except Exception as e:
        print(f"Error in report generation: {e}")
        return False


def process_vcf_complete(vcf_path, output_dir, sample_id="sample", analysis_result=None):
    """Complete pipeline: VCF -> CSV -> Features -> Scores -> Report"""
    
    print(f"DEBUG: Starting VCF processing for {vcf_path}")
    print(f"DEBUG: Output directory: {output_dir}")
    
    # If VCF content is in database, write it to temp file for processing
    if analysis_result and analysis_result.upload.vcf_content:
        print(f"DEBUG: Using VCF content from database")
        temp_vcf_path = os.path.join(output_dir, f"{sample_id}.vcf")
        # Ensure output directory exists before writing temp file
        os.makedirs(output_dir, exist_ok=True)
        with open(temp_vcf_path, 'w', encoding='utf-8') as f:
            f.write(analysis_result.upload.vcf_content)
        vcf_path = temp_vcf_path
        print(f"DEBUG: Wrote VCF content to temporary file: {vcf_path}")
    else:
        print(f"DEBUG: Using VCF file from filesystem")
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Check if input file exists
    if not os.path.exists(vcf_path):
        return None, f"Input file does not exist: {vcf_path}"
    
    print(f"DEBUG: Input file exists, size: {os.path.getsize(vcf_path)} bytes")
    
    # Handle VCF.gz files - decompress first
    if vcf_path.endswith('.gz'):
        import gzip
        decompressed_path = vcf_path[:-3]  # Remove .gz extension
        try:
            with gzip.open(vcf_path, 'rt') as gz_file:
                with open(decompressed_path, 'w') as vcf_file:
                    vcf_file.write(gz_file.read())
            print(f"DEBUG: Decompressed VCF.gz to {decompressed_path}")
            vcf_path = decompressed_path
        except Exception as e:
            print(f"DEBUG: Decompression failed: {e}")
            return None, f"Failed to decompress VCF.gz: {e}"
    
    # Step 1: VCF to CSV
    csv_path = os.path.join(output_dir, "variants.csv")
    print(f"DEBUG: Converting VCF to CSV: {vcf_path} -> {csv_path}")
    
    if not parse_vcf_to_csv(vcf_path, csv_path):
        print(f"DEBUG: VCF to CSV conversion failed")
        return None, "Failed to convert VCF to CSV"
    
    print(f"DEBUG: VCF to CSV conversion successful, CSV size: {os.path.getsize(csv_path)} bytes")
    
    # Step 2: CSV to Features
    features_path = os.path.join(output_dir, "features.csv")
    if not process_variant_features(csv_path, features_path):
        return None, "Failed to extract features"
    
    # Step 3: Features to Scores
    scores_path = os.path.join(output_dir, "scores.csv")
    if not apply_ai_scoring(features_path, scores_path):
        return None, "Failed to apply AI scoring"
    
    # Step 4: Scores to Report
    report_path = os.path.join(output_dir, "clinical_report.html")
    if not generate_report(scores_path, report_path, sample_id):
        return None, "Failed to generate report"
    
    # Step 5: Populate database with generated content
    if analysis_result:
        try:
            # Read scores CSV and store safely as JSON compatible structures
            import json
            scores_df = pd.read_csv(scores_path)
            
            # to_json handles all NaN/inf float formatting natively, json.loads produces pristine None items
            clean_json_str = scores_df.to_json(orient='records')
            analysis_result.scores_data = json.loads(clean_json_str)
            
            # Read HTML report
            with open(report_path, 'r', encoding='utf-8') as f:
                analysis_result.clinical_report_html = f.read()
            
            # Read and encode charts as base64
            import base64
            
            chart_paths = {
                'gene_roles_chart': report_path.replace('.html', '_gene_roles.png'),
                'mutation_origin_chart': report_path.replace('.html', '_mutation_origin.png'),
                'vaf_distribution_chart': report_path.replace('.html', '_vaf_dist.png')
            }
            
            for field_name, chart_path in chart_paths.items():
                if os.path.exists(chart_path):
                    with open(chart_path, 'rb') as f:
                        chart_data = base64.b64encode(f.read()).decode('utf-8')
                        setattr(analysis_result, field_name, chart_data)
            
            analysis_result.analysis_status = "completed"
            analysis_result.save()
            
            print(f"DEBUG: Successfully populated database with analysis results")
            
            # Clean up temporary files to save disk space in production
            try:
                import shutil
                if os.path.exists(output_dir):
                    shutil.rmtree(output_dir)
                    print(f"DEBUG: Cleaned up temporary files from {output_dir}")
            except Exception as cleanup_error:
                print(f"DEBUG: Cleanup warning: {cleanup_error}")
                # Don't fail the whole process if cleanup fails
            
        except Exception as e:
            print(f"DEBUG: Failed to populate database: {e}")
            return None, f"Failed to save results to database: {e}"
    
    return report_path, "Success"
