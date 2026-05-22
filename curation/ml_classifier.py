# =============================================================================
# CLASSIFICATEUR IA POUR VARIANTES GÉNOMIQUES
# =============================================================================
# Ce module utilise un Random Forest pour classifier les variants génomiques
# en fonction de leur probabilité d'être pathogènes.
#
# Fonctionnalités principales:
# - Extraction de caractéristiques depuis les données VCF
# - Encodage des rôles de gènes et origines de mutations
# - Classification avec Random Forest (100 arbres)
# - Scoring de probabilité pathogène (0-1)

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
import pickle
import os

class SimpleAIClassifier:
    """Classificateur IA simple utilisant Random Forest pour les variants génomiques"""
    
    def __init__(self):
        # Initialiser le modèle Random Forest avec 100 arbres
        self.model = RandomForestClassifier(n_estimators=100, random_state=42)
        
        # Encodeurs pour les caractéristiques catégorielles
        self.gene_role_encoder = LabelEncoder()  # Rôle du gène (Oncogène, etc.)
        self.origin_encoder = LabelEncoder()   # Origine de la mutation (Somatic, Germline)
        self.is_trained = False
        
    def extract_basic_features(self, df):
        """Extraire les caractéristiques de base depuis les données de variants"""
        features = []
        
        for _, row in df.iterrows():
            # Caractéristiques numériques de base
            vaf = float(row.get('vaf', 0))  # Variant Allele Frequency
            depth = float(row.get('depth', 1))  # Profondeur de séquençage
            log_depth = np.log10(depth)  # Échelle log pour la profondeur
            
            # Encodage du rôle du gène
            gene_role = str(row.get('gene_role', 'Unknown'))
            try:
                gene_role_encoded = self.gene_role_encoder.transform([gene_role])[0]
            except ValueError:
                gene_role_encoded = self.gene_role_encoder.transform(['Unknown'])[0]
            
            # Encodage de l'origine de la mutation
            origin = str(row.get('origin_improved', row.get('origin', 'Unknown')))
            try:
                origin_encoded = self.origin_encoder.transform([origin])[0]
            except ValueError:
                origin_encoded = self.origin_encoder.transform(['Unknown'])[0]
            
            # Statut LOH (Loss of Heterozygosity)
            loh_status = 1 if 'LOH' in str(row.get('loh_status', '')) else 0
            
            # Combiner toutes les caractéristiques
            feature_vector = [vaf, depth, log_depth, gene_role_encoded, origin_encoded, loh_status]
            features.append(feature_vector)
            
        return np.array(features)
    
    def train(self, df):
        """Entraîner le modèle avec des étiquettes synthétiques simples"""
        print("Entraînement du modèle IA simple...")
        
        # Préparer les encodeurs
        gene_roles = list(df['gene_role'].unique()) + ['Unknown']
        if 'origin' in df.columns:
            origins = list(df['origin_improved'].fillna(df['origin']).unique()) + ['Unknown']
        else:
            origins = list(df['origin_improved'].fillna('Unknown').unique()) + ['Unknown']
        
        self.gene_role_encoder.fit(gene_roles)
        self.origin_encoder.fit(origins)
        
        # Extraire les caractéristiques
        X = self.extract_basic_features(df)
        
        # Créer des étiquettes synthétiques équilibrées basées sur les caractéristiques des variants
        scores = []
        for _, row in df.iterrows():
            score = 0
            vaf = row.get('vaf', 0)
            depth = row.get('depth', 0)
            gene_role = str(row.get('gene_role', ''))
            origin = str(row.get('origin_improved', row.get('origin', '')))
            mut_status = row.get('mut_status', '')
            gene_val = str(row.get('gene', '')).strip()
            
            # Règles de classification définies par l'utilisateur (priorité la plus élevée)
            # Règle 1: Intergénique -> Signification clinique inconnue (priorité la plus basse)
            if gene_val == 'Intergenic':
                score = 0.1  # Très faible priorité
                scores.append(score)
                continue
                
            # Règle 2: Variante somatique confirmée -> Signification clinique très forte (priorité la plus élevée)
            if pd.notna(mut_status) and 'Confirmed somatic variant' in str(mut_status):
                score = 1.0  # Très haute priorité
                scores.append(score)
                continue
                
            # Règle 3: Statut NaN -> Signification clinique potentielle (priorité moyenne)
            if pd.isna(mut_status) or mut_status == '':
                score = 0.5  # Priorité moyenne
                scores.append(score)
                continue
            
            # Scoring traditionnel pour les variants restants
            # Scoring VAF (contribution continue)
            if vaf > 0.5:
                score += 0.4
            elif vaf > 0.2:
                score += 0.3
            elif vaf > 0.05:
                score += 0.2
            else:
                score += 0.1
                
            # Depth scoring
            if depth > 200:
                score += 0.2
            elif depth > 50:
                score += 0.15
            else:
                score += 0.05
                
            # Gene role scoring
            if 'Oncogene' in gene_role:
                score += 0.25
            elif 'Tumor Suppressor' in gene_role or 'TSG' in gene_role:
                score += 0.3
            else:
                score += 0.05
                
            # Origin scoring
            if 'Somatic' in origin and 'LOH' in origin:
                score += 0.15
            elif 'Somatic' in origin:
                score += 0.1
            elif 'Germline' in origin:
                score += 0.05
            else:
                score += 0.02
                
            scores.append(score)
        
        # Create binary labels: 1 for high priority (score > 0.5), 0 for low
        y = [1 if s > 0.5 else 0 for s in scores]
        
        # Train model
        self.model.fit(X, y)
        self.is_trained = True
        
        print(f"Model trained on {len(X)} variants")
        print(f"High priority variants: {sum(y)} out of {len(y)} ({sum(y)/len(y)*100:.1f}%)")
        
    def predict(self, df):
        """Make predictions with user-defined classification rules"""
        if not self.is_trained:
            self.train(df)
        
        X = self.extract_basic_features(df)
        probabilities = self.model.predict_proba(X)  # Get probabilities for both classes
        ai_predictions = probabilities[:, 1]  # Probability of class 1 (high priority)
        
        result_df = df.copy()
        result_df['ai_score'] = ai_predictions
        
        # Apply user-defined classification logic
        def classify_variant(row):
            mut_status = row.get('mut_status', '')
            gene_val = str(row.get('gene', '')).strip()
            
            # Rule 1: Intergenic -> Signification clinique inconnue (highest precedence)
            if gene_val == 'Intergenic':
                return 0.1  # Lowest priority
            
            # Rule 2: Confirmed somatic variant -> Signification clinique très forte
            elif pd.notna(mut_status) and 'Confirmed somatic variant' in str(mut_status):
                return 1.0  # Highest priority
            
            # Rule 3: Status NaN -> Signification clinique potentielle
            elif pd.isna(mut_status) or mut_status == '':
                return 0.5  # Medium priority
            
            # Fallback: AI-based priority classification
            else:
                return row.get('ai_score', 0)
        
        # Apply classification and update priority_score
        result_df['priority_score'] = result_df.apply(classify_variant, axis=1)
        
        return result_df
    
    def save_model(self, filepath):
        """Save the trained model"""
        model_data = {
            'model': self.model,
            'gene_role_encoder': self.gene_role_encoder,
            'origin_encoder': self.origin_encoder,
            'is_trained': self.is_trained
        }
        
        with open(filepath, 'wb') as f:
            pickle.dump(model_data, f)
        print(f"Model saved to {filepath}")
    
    def load_model(self, filepath):
        """Load a trained model"""
        with open(filepath, 'rb') as f:
            model_data = pickle.load(f)
        
        self.model = model_data['model']
        self.gene_role_encoder = model_data['gene_role_encoder']
        self.origin_encoder = model_data['origin_encoder']
        self.is_trained = model_data['is_trained']
        
        print(f"Model loaded from {filepath}")


def main():
    """Main execution"""
    try:
        # Load input data
        df = pd.read_csv(snakemake.input.features)
        
        if df.empty:
            pd.DataFrame().to_csv(snakemake.output.scores, index=False)
            return
        
        # Initialize simple AI classifier
        classifier = SimpleAIClassifier()
        
        # Check for existing model
        model_path = "/home/zakaria/variant_curation_app/models/simple_ai_model.pkl"
        
        if os.path.exists(model_path):
            print("Loading existing AI model...")
            classifier.load_model(model_path)
            result_df = classifier.predict(df)
        else:
            print("Training new AI model...")
            classifier.train(df)
            result_df = classifier.predict(df)
            
            # Save model
            os.makedirs("/home/zakaria/variant_curation_app/models", exist_ok=True)
            classifier.save_model(model_path)
        
        # Save results
        result_df.to_csv(snakemake.output.scores, index=False)
        
        print(f"AI processing complete for {len(result_df)} variants")
        print(f"Mean priority score: {result_df['priority_score'].mean():.3f}")
        print(f"Priority score range: {result_df['priority_score'].min():.3f} - {result_df['priority_score'].max():.3f}")
        print(f"High confidence variants: {(result_df['priority_score'] > 0.7).sum()}")
        
    except Exception as e:
        print(f"AI Error: {e}")
        raise e


if __name__ == "__main__":
    main()
