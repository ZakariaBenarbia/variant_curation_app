import zipfile
import io
import os
from django.http import HttpResponse
from django.shortcuts import get_object_or_404
from .models import VCFUpload, AnalysisResult
import json
import base64
import pandas as pd

def download_results_zip(request, upload_id):
    """Download all results from database for an upload as a ZIP file"""
    upload = get_object_or_404(VCFUpload, pk=upload_id, user=request.user)
    
    # Create a ZIP file in memory
    zip_buffer = io.BytesIO()
    
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        # Add clinical report HTML
        if upload.analysis_result and upload.analysis_result.clinical_report_html:
            zip_file.writestr(f"{upload.title}_clinical_report.html", 
                           upload.analysis_result.clinical_report_html)
        
        # Add scores data as CSV
        if upload.analysis_result and upload.analysis_result.scores_data:
            df = pd.DataFrame(upload.analysis_result.scores_data)
            csv_buffer = io.StringIO()
            df.to_csv(csv_buffer, index=False)
            zip_file.writestr(f"{upload.title}_scores.csv", csv_buffer.getvalue())
        
        # Add charts as PNG files
        if upload.analysis_result:
            charts = [
                ('gene_roles_chart', upload.analysis_result.gene_roles_chart),
                ('mutation_origin_chart', upload.analysis_result.mutation_origin_chart), 
                ('vaf_distribution_chart', upload.analysis_result.vaf_distribution_chart)
            ]
            
            for chart_name, chart_data in charts:
                if chart_data:
                    try:
                        # Decode base64 and add as PNG
                        image_data = base64.b64decode(chart_data)
                        zip_file.writestr(f"{upload.title}_{chart_name}.png", image_data)
                    except Exception as e:
                        print(f"Error decoding {chart_name}: {e}")
        
        # Add original VCF file content
        if upload.vcf_content:
            zip_file.writestr(f"{upload.title}.vcf", upload.vcf_content)
        
        # Add metadata
        metadata = {
            'upload_id': upload.id,
            'title': upload.title,
            'uploaded_at': str(upload.uploaded_at),
            'sample_name': upload.sample_name,
            'analysis_status': upload.analysis_result.analysis_status if upload.analysis_result else 'N/A',
            'created_at': str(upload.analysis_result.created_at) if upload.analysis_result else 'N/A',
            'updated_at': str(upload.analysis_result.updated_at) if upload.analysis_result else 'N/A'
        }
        
        zip_file.writestr(f"{upload.title}_metadata.json", json.dumps(metadata, indent=2))
    
    # Prepare response
    zip_buffer.seek(0)
    response = HttpResponse(zip_buffer, content_type='application/zip')
    response['Content-Disposition'] = f'attachment; filename="{upload.title}_results.zip"'
    
    return response
