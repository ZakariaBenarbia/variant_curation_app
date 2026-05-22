from django.db import models
from django.contrib.auth.models import User
from django.utils.text import slugify
from django.db.models.signals import post_save
from django.dispatch import receiver

class TeamMember(models.Model):
    name = models.CharField(max_length=100)
    email = models.EmailField()
    phone = models.CharField(max_length=20, blank=True)
    # Password to be set by admin when creating the team member
    password = models.CharField(max_length=128, blank=True, help_text="Set the login password for this team member")
    # Link to Django's User model for secure authentication
    user = models.OneToOneField(User, null=True, blank=True, on_delete=models.CASCADE, related_name='team_member')

    def __str__(self):
        return self.name

class VCFUpload(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    title = models.CharField(max_length=200)
    vcf_file = models.FileField(upload_to="vcf_uploads/", blank=True, null=True)
    vcf_content = models.TextField(blank=True, null=True, help_text="VCF file content stored in database")
    uploaded_at = models.DateTimeField(auto_now_add=True)
    sample_name = models.CharField(max_length=200, blank=True)

    def __str__(self):
        return f"{self.title} ({self.uploaded_at:%Y-%m-%d %H:%M})"


class AnalysisResult(models.Model):
    """Store analysis results that were previously saved as files"""
    upload = models.OneToOneField(VCFUpload, on_delete=models.CASCADE, related_name="analysis_result")
    scores_data = models.JSONField(default=dict, help_text="Complete scores CSV data")
    clinical_report_html = models.TextField(blank=True, help_text="Generated clinical report HTML")
    gene_roles_chart = models.TextField(blank=True, help_text="Base64 encoded gene roles chart")
    mutation_origin_chart = models.TextField(blank=True, help_text="Base64 encoded mutation origin chart")
    vaf_distribution_chart = models.TextField(blank=True, help_text="Base64 encoded VAF distribution chart")
    analysis_status = models.CharField(max_length=50, default="pending")
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self):
        return f"Analysis for {self.upload.title}"

class ProcessingLog(models.Model):
    """Track processing steps and logs"""
    upload = models.ForeignKey(VCFUpload, on_delete=models.CASCADE, related_name="processing_logs")
    step = models.CharField(max_length=100)
    status = models.CharField(max_length=50)  # success, error, warning
    message = models.TextField()
    timestamp = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ["-timestamp"]

    def __str__(self):
        return f"{self.upload.title} - {self.step} ({self.status})"

# Signal: Create or update a corresponding Django User when TeamMember is saved
@receiver(post_save, sender=TeamMember)
def create_or_update_user(sender, instance, created, **kwargs):
    """
    After a TeamMember is saved, ensure a corresponding User exists.
    Username = slugified name (unique). User is linked back to TeamMember.
    """
    # If user already linked, nothing to do
    if instance.user:
        return

    # Generate unique username from team member name
    base_username = slugify(instance.name) if instance.name else f"user_{instance.id}"
    username = base_username
    counter = 1

    # Ensure uniqueness by appending counter if needed
    while User.objects.filter(username=username).exists():
        username = f"{base_username}{counter}"
        counter += 1

    # Create User with the generated username
    user = User.objects.create_user(
        username=username,
        email=instance.email if instance.email else '',
        first_name=instance.name
    )
    # Password will be set by admin via the password field in TeamMember
    if instance.password:
        user.set_password(instance.password)
        user.save()
    user.save()
    # Link the User back to the TeamMember
    instance.user = user
    instance.save(update_fields=['user'])
