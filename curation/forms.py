from django import forms
from django.contrib.auth.forms import AuthenticationForm

class CustomAuthenticationForm(AuthenticationForm):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['username'].widget.attrs.update({'placeholder': "Nom d'utilisateur"})
        self.fields['password'].widget.attrs.update({'placeholder': 'Mot de passe'})
    
    error_messages = {
        'invalid_login': "Identifiant ou mot de passe incorrect.",
        'inactive': "Ce compte est inactif.",
        'invalid': "Veuillez entrer un identifiant et un mot de passe valides.",
        'required': "Ce champ est obligatoire.",
    }

class VCFUploadForm(forms.Form):
    vcf_file = forms.FileField(label="Fichier VCF", help_text="Importer un fichier .vcf ou .vcf.gz")
