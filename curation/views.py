import os
import shutil
import pandas as pd
import base64
from django.conf import settings
from django.http import FileResponse, Http404, HttpResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse
from django.core.files.storage import FileSystemStorage
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import User
from .forms import VCFUploadForm
from .models import VCFUpload, AnalysisResult, ProcessingLog, TeamMember
from .django_adapters import process_vcf_complete
from .download_results import download_results_zip


def _delete_upload_data(upload):
    """Delete upload files and database records."""
    # Delete analysis results from database
    if hasattr(upload, 'analysis_result'):
        upload.analysis_result.delete()
    
    # Delete processing logs
    upload.processing_logs.all().delete()

    # Clear VCF content from database
    upload.vcf_content = None
    upload.save(update_fields=['vcf_content'])

    processed_dir = os.path.join(settings.MEDIA_ROOT, f"processed_{upload.id}")
    if os.path.exists(processed_dir):
        shutil.rmtree(processed_dir)

    if upload.vcf_file and os.path.exists(upload.vcf_file.path):
        os.remove(upload.vcf_file.path)

    upload.delete()


@login_required
def upload_vcf(request):
    uploads = VCFUpload.objects.filter(user=request.user).order_by("-uploaded_at")[:10]
    if request.method == "POST":
        form = VCFUploadForm(request.POST, request.FILES)
        if form.is_valid():
            vcf_file = request.FILES["vcf_file"]
            original_name = os.path.basename(vcf_file.name)
            if original_name.lower().endswith(".vcf.gz"):
                title = original_name[:-7]
            elif original_name.lower().endswith(".vcf"):
                title = original_name[:-4]
            else:
                title = os.path.splitext(original_name)[0]
            
            # Read VCF file content for database storage
            if original_name.lower().endswith(".vcf.gz"):
                # Handle compressed VCF files
                import gzip
                with gzip.open(vcf_file, 'rt', encoding='utf-8') as gz_file:
                    vcf_content = gz_file.read()
            else:
                # Handle regular VCF files
                vcf_content = vcf_file.read().decode('utf-8')
            
            # Save VCF file temporarily for processing
            fs = FileSystemStorage(location=settings.MEDIA_ROOT)
            filename = fs.save(vcf_file.name, vcf_file)
            file_path = fs.path(filename)
            
            # Create upload record with VCF content in database
            upload = VCFUpload.objects.create(
                user=request.user,
                vcf_file=filename,
                title=title,
                vcf_content=vcf_content,
                sample_name=title,
            )
            
            # Create analysis result record
            analysis_result = AnalysisResult.objects.create(upload=upload)
            ProcessingLog.objects.create(upload=upload, step="upload", status="success", message="VCF file uploaded successfully")
            
            # Use the new Django adapters to process VCF file
            output_dir = os.path.join(settings.MEDIA_ROOT, f"processed_{upload.id}")
            report_path, status = process_vcf_complete(file_path, output_dir, upload.title, analysis_result)
            
            if status != "Success":
                analysis_result.analysis_status = "failed"
                analysis_result.save()
                ProcessingLog.objects.create(upload=upload, step="processing", status="error", message=f"Processing failed: {status}")
                return render(request, "curation/upload.html", {
                    "form": form,
                    "uploads": VCFUpload.objects.filter(user=request.user).order_by("-uploaded_at")[:10],
                    "error": f"Processing failed: {status}"
                })
            
            ProcessingLog.objects.create(upload=upload, step="processing", status="success", message="VCF processing completed")
                
            return redirect(reverse("curation:upload_vcf"))
    else:
        form = VCFUploadForm()
    return render(request, "curation/upload.html", {"form": form, "uploads": uploads})




@login_required
def download_upload_csv(request, upload_id):
    """Download processed CSV data from database for an upload"""
    upload = get_object_or_404(VCFUpload, pk=upload_id, user=request.user)
    
    # Try to get analysis result from database
    if hasattr(upload, 'analysis_result') and upload.analysis_result.scores_data:
        # Convert JSON data back to CSV
        df = pd.DataFrame(upload.analysis_result.scores_data)
        
        # Create CSV response
        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = f'attachment; filename="{upload.title}_scores.csv"'
        
        # Write DataFrame to CSV
        df.to_csv(response, index=False)
        return response
    else:
        return HttpResponse("Analysis data not found in database", status=404)


@login_required
def upload_detail(request, upload_id):
    """Display upload details - now redirects to upload page since variants are stored in analysis_result"""
    upload = get_object_or_404(VCFUpload, pk=upload_id, user=request.user)
    
    # Template removed - redirect to upload page instead
    return redirect(reverse("curation:upload_vcf"))



@login_required
def delete_upload(request, upload_id):
    """Delete an upload and all associated data"""
    upload = get_object_or_404(VCFUpload, pk=upload_id, user=request.user)
    
    if request.method == "POST":
        _delete_upload_data(upload)
        return redirect(reverse("curation:upload_vcf"))
    
    return render(request, "curation/delete_confirm.html", {"upload": upload})


@login_required
def bulk_delete_uploads(request):
    """Delete selected uploads for the current user."""
    if request.method != "POST":
        return redirect(reverse("curation:upload_vcf"))

    upload_ids = request.POST.getlist("upload_ids")
    if upload_ids:
        uploads = VCFUpload.objects.filter(user=request.user, id__in=upload_ids)
        for upload in uploads:
            _delete_upload_data(upload)

    return redirect(reverse("curation:upload_vcf"))


@login_required
def download_results_zip_view(request, upload_id):
    """Download all results for an upload as a ZIP file"""
    return download_results_zip(request, upload_id)


@login_required
def generate_report(request, upload_id):
    """Serve clinical report from database for an upload"""
    upload = get_object_or_404(VCFUpload, pk=upload_id, user=request.user)
    
    # Try to get analysis result from database
    if hasattr(upload, 'analysis_result') and upload.analysis_result.clinical_report_html:
        content = upload.analysis_result.clinical_report_html
        response = HttpResponse(content, content_type="text/html")
        return response
    else:
        return HttpResponse("Clinical report not found in database. Please process the upload first.", status=404)


@login_required
def report_asset(request, upload_id, asset_name):
    """Serve generated report image assets from database for an upload."""
    allowed_assets = {
        "clinical_report_gene_roles.png": "gene_roles_chart",
        "clinical_report_mutation_origin.png": "mutation_origin_chart", 
        "clinical_report_vaf_dist.png": "vaf_distribution_chart",
    }
    
    if asset_name not in allowed_assets:
        raise Http404("Asset not found")

    upload = get_object_or_404(VCFUpload, pk=upload_id, user=request.user)
    
    # Get chart data from database
    if hasattr(upload, 'analysis_result'):
        field_name = allowed_assets[asset_name]
        chart_data = getattr(upload.analysis_result, field_name, None)
        
        if chart_data:
            # Decode base64 and return as image
            try:
                image_data = base64.b64decode(chart_data)
                return HttpResponse(image_data, content_type="image/png")
            except Exception as e:
                return HttpResponse(f"Error decoding chart: {e}", status=500)
    
    raise Http404("Report asset not found in database. Generate report first.")


@login_required
def add_team_member(request):
    """Add new team member view for admins."""
    if not request.user.is_staff:
        return redirect('upload_vcf')
    
    if request.method == 'POST':
        name = request.POST.get('name', '').strip()
        email = request.POST.get('email', '').strip()
        phone = request.POST.get('phone', '').strip()
        password = request.POST.get('password', '').strip()
        
        if name and email:
            # Create team member - the signal will automatically create the Django User
            team_member = TeamMember.objects.create(
                name=name,
                email=email,
                phone=phone,
                password=password
            )
            
            return render(request, 'curation/team_member_added.html', {'member': team_member})
    
    return render(request, 'curation/add_team_member.html')


@login_required
def team_members_list(request):
    """Show all team members for admin management."""
    if not request.user.is_staff:
        return redirect('upload_vcf')
    
    team_members = TeamMember.objects.all().order_by('name')
    return render(request, 'curation/team_members_list.html', {'team_members': team_members})


@login_required
def delete_team_member(request, member_id):
    """Delete team member for admins."""
    if not request.user.is_staff:
        return redirect('curation:team_members_list')
    
    try:
        member = TeamMember.objects.get(id=member_id)
        
        if request.method == 'POST':
            member.delete()
            return redirect('curation:team_members_list')
        else:
            return render(request, 'curation/team_member_delete_confirm.html', {'member': member})
            
    except TeamMember.DoesNotExist:
        return redirect('curation:team_members_list')
