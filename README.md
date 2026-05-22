# Variant Curation Django App

This folder contains a lightweight Django app scaffold for generic VCF ingestion and review.

## Setup

1. Create and activate a Python environment.
2. Install dependencies:

```bash
pip install -r requirements.txt
```

3. Run migrations:

```bash
python manage.py migrate
```

4. Start the server:

```bash
python manage.py runserver
```

5. Open `http://127.0.0.1:8000/` to upload a VCF and review parsed variants.

## Notes

- The app currently parses any VCF via `pysam` and imports standard fields.
- `Variant` records store generic VCF values plus `vaf`, `depth`, and `allelic_depth` when available.
- This is a starting scaffold for variant curation and review.
