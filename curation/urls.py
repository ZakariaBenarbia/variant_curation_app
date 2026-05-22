from django.urls import path
from . import views

app_name = "curation"

urlpatterns = [
    path("", views.upload_vcf, name="upload_vcf"),
    path("upload/<int:upload_id>/", views.upload_detail, name="upload_detail"),
    path("upload/<int:upload_id>/delete/", views.delete_upload, name="delete_upload"),
    path("uploads/bulk-delete/", views.bulk_delete_uploads, name="bulk_delete_uploads"),
    path("upload/<int:upload_id>/download/", views.download_upload_csv, name="download_upload_csv"),
    path("upload/<int:upload_id>/download-results/", views.download_results_zip_view, name="download_results"),
    path("upload/<int:upload_id>/report/", views.generate_report, name="generate_report"),
    path("upload/<int:upload_id>/report-asset/<str:asset_name>/", views.report_asset, name="report_asset"),
    path("add-team-member/", views.add_team_member, name="add_team_member"),
    path("team-members/", views.team_members_list, name="team_members_list"),
    path("team-members/<int:member_id>/delete/", views.delete_team_member, name="delete_team_member"),
]
