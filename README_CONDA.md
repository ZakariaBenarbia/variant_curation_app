# PostgreSQL Setup with Conda

## 1. Create Conda Environment
```bash
conda env create -f environment.yml
conda activate variant-curation-app
```

## 2. Install PostgreSQL (if not already installed)

### Option A: Install via system package manager (recommended)
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install postgresql postgresql-contrib

# macOS with Homebrew
brew install postgresql
brew services start postgresql

# Windows
# Download from https://www.postgresql.org/download/windows/
```

### Option B: Install via conda (for development only)
```bash
conda install -c conda-forge postgresql
```

## 3. Set Up Database
```bash
# Start PostgreSQL service
sudo service postgresql start  # Linux
brew services start postgresql  # macOS

# Create database and user
sudo -u postgres psql
CREATE DATABASE variant_curation_db;
CREATE USER your_username WITH PASSWORD 'your_password';
GRANT ALL PRIVILEGES ON DATABASE variant_curation_db TO your_username;
\q
```

## 4. Configure Environment
```bash
# Copy environment template
cp .env.example .env

# Edit .env file with your database credentials
nano .env
```

## 5. Run Django Migrations
```bash
python manage.py makemigrations
python manage.py migrate
python manage.py createsuperuser
```

## 6. Test the Setup
```bash
python manage.py runserver
```

## Notes
- The conda environment includes psycopg2 (PostgreSQL adapter)
- python-decouple is installed via pip since it's not available in conda-forge
- pysam is installed via pip as it's system-dependent
