from django.contrib import admin
from .models import VCFUpload, TeamMember, AnalysisResult, ProcessingLog

@admin.register(TeamMember)
class TeamMemberAdmin(admin.ModelAdmin):
    list_display = ("name", "email", "phone", "user")
    search_fields = ("name", "email", "phone")
    ordering = ("name",)
    readonly_fields = ("user",)  # User is created automatically; prevent manual changes here
    fields = ("name", "email", "phone", "password", "user")  # Show password field in form

@admin.register(VCFUpload)
class VCFUploadAdmin(admin.ModelAdmin):
    list_display = ("title", "user", "sample_name", "uploaded_at")
    list_filter = ("user", "uploaded_at")
    search_fields = ("title", "sample_name", "user__username")
    ordering = ("-uploaded_at",)

@admin.register(AnalysisResult)
class AnalysisResultAdmin(admin.ModelAdmin):
    list_display = ("upload", "analysis_status", "created_at", "updated_at")
    list_filter = ("analysis_status", "created_at")
    search_fields = ("upload__title", "upload__user__username")
    ordering = ("-created_at",)

@admin.register(ProcessingLog)
class ProcessingLogAdmin(admin.ModelAdmin):
    list_display = ("upload", "step", "status", "timestamp")
    list_filter = ("status", "step", "timestamp")
    search_fields = ("upload__title", "step", "message")
    ordering = ("-timestamp",)
