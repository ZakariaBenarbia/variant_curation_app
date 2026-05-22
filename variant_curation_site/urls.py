from django.contrib import admin
from django.urls import include, path
from django.conf import settings
from django.conf.urls.static import static
from django.contrib.auth import views as auth_views
from django.views.generic import RedirectView
from curation.forms import CustomAuthenticationForm

urlpatterns = [
    path("admin/", admin.site.urls),
    # Redirect root URL to curation app
    path("", RedirectView.as_view(pattern_name="curation:upload_vcf", permanent=False), name="home"),
    # Login Route
    path("login/", auth_views.LoginView.as_view(template_name="curation/login.html", authentication_form=CustomAuthenticationForm), name="login"),
    # Logout Route: redirect to login page after logout
    path("logout/", auth_views.LogoutView.as_view(next_page="/login/"), name="logout"),
    path("curation/", include("curation.urls", namespace="curation")),
]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
