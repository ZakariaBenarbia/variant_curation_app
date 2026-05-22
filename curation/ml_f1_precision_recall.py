#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
MÉTRIQUES DE CLASSIFICATION : F1-SCORE, PRECISION, RECALL
================================================================================
Évaluation simple du modèle Random Forest avec métriques standard.

Usage :
    python ml_f1_precision_recall.py
    python ml_f1_precision_recall.py <fichier_scores.csv>

Auteur : Zakaria
Date : Avril 2026
================================================================================
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_curve, auc
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
import os
import sys


def generate_binary_labels(df):
    """
    Génère les labels binaires (0/1) basés sur les règles expertes
    0: Basse priorité (score <= 0.5)
    1: Haute priorité (score > 0.5)
    """
    labels = []
    
    for _, row in df.iterrows():
        gene_val = str(row.get('gene', '')).strip()
        mut_status = str(row.get('mut_status', ''))
        
        # Classe 1: Haute priorité (COSMIC confirmé)
        if 'Confirmed somatic variant' in mut_status:
            labels.append(1)
        # Classe 0: Basse priorité (Intergenic)
        elif gene_val == 'Intergenic':
            labels.append(0)
        # Classe 1: Haute priorité (mut_status vide)
        elif mut_status == '' or mut_status == 'nan':
            labels.append(1)
        else:
            # Pour les autres, utiliser improved_score comme heuristique
            score = float(row.get('improved_score', 0))
            if score > 0.5:
                labels.append(1)
            else:
                labels.append(0)
    
    return np.array(labels)


def prepare_features(df):
    """Prépare les features biologiques (6 features)"""
    gene_role_encoder = LabelEncoder()
    origin_encoder = LabelEncoder()
    
    gene_roles = list(df['gene_role'].unique()) + ['Unknown']
    origins = list(df['origin_improved'].fillna('Unknown').unique()) + ['Unknown']
    
    gene_role_encoder.fit(gene_roles)
    origin_encoder.fit(origins)
    
    X = []
    for _, row in df.iterrows():
        # Features numériques de base
        vaf = float(row.get('vaf', 0))
        depth = float(row.get('depth', 1))
        log_depth = np.log10(depth) if depth > 0 else 0
        
        # Encodage catégoriel
        gene_role = str(row.get('gene_role', 'Unknown'))
        try:
            gene_role_enc = gene_role_encoder.transform([gene_role])[0]
        except ValueError:
            gene_role_enc = gene_role_encoder.transform(['Unknown'])[0]
        
        origin = str(row.get('origin_improved', 'Unknown'))
        try:
            origin_enc = origin_encoder.transform([origin])[0]
        except ValueError:
            origin_enc = origin_encoder.transform(['Unknown'])[0]
        
        loh = 1 if 'LOH' in str(row.get('loh_status', '')) else 0
        
        # 6 features biologiques
        X.append([vaf, depth, log_depth, gene_role_enc, origin_enc, loh])
    
    return np.array(X)


def main():
    """Point d'entrée principal"""
    
    # Load the dataset
    if len(sys.argv) > 1:
        data_path = sys.argv[1]
    else:
        candidates = [
            "scores.csv",
            "/home/zakaria/Bureau/Nouveau dossier/SRR29289914_filtered_scores.csv"
        ]
        data_path = None
        for c in candidates:
            if os.path.exists(c):
                data_path = c
                break
        
        if not data_path:
            print("Usage: python ml_f1_precision_recall.py <fichier_scores.csv>")
            sys.exit(1)
    
    print(f"Chargement du dataset : {data_path}")
    variant_data = pd.read_csv(data_path)
    print(f"✓ {len(variant_data)} variants chargés\n")
    
    # Générer les labels binaires basés sur les règles expertes
    y = generate_binary_labels(variant_data)
    
    # Préparer les features biologiques
    X = prepare_features(variant_data)
    
    print("=" * 60)
    print("ÉVALUATION DU MODÈLE ML BINAIRE (0/1)")
    print("=" * 60)
    print(f"\nFeatures biologiques : {X.shape[1]} dimensions")
    print(f"  - VAF, depth, log_depth")
    print(f"  - gene_role_encoded, origin_encoded, LOH")
    print(f"\nLabels binaires basés sur les règles expertes:")
    print(f"  - COSMIC confirmé → classe 1 (haute priorité)")
    print(f"  - Intergenic → classe 0 (basse priorité)")
    print(f"  - mut_status vide → classe 1 (haute priorité)")
    print(f"  - improved_score > 0.5 → classe 1")
    print(f"  - improved_score <= 0.5 → classe 0")
    print(f"\nDistribution des classes:")
    print(f"  - Basse priorité (0): {sum(y==0)} ({sum(y==0)/len(y)*100:.1f}%)")
    print(f"  - Haute priorité (1): {sum(y==1)} ({sum(y==1)/len(y)*100:.1f}%)")
    print("=" * 60 + "\n")
    
    # Split the dataset into training (70%) and testing (30%) sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42, stratify=y)
    
    print(f"Train : {len(X_train)} échantillons")
    print(f"Test  : {len(X_test)} échantillons\n")
    
    # Initialize the Random Forest classifier (binary)
    random_forest_model = RandomForestClassifier(n_estimators=100, random_state=42)
    
    # Train the model on the training data
    print("Entraînement du modèle...")
    random_forest_model.fit(X_train, y_train)
    print("✓ Entraînement terminé\n")
    
    # Predictions on the test data
    y_pred = random_forest_model.predict(X_test)
    
    # Evaluate the model's performance
    accuracy = accuracy_score(y_test, y_pred)
    classification_rep = classification_report(y_test, y_pred, 
                                               target_names=['Basse priorité', 'Haute priorité'],
                                               digits=3)
    
    print("=" * 60)
    print("RÉSULTATS DE L'ÉVALUATION")
    print("=" * 60)
    print(f'\nAccuracy = {accuracy:.4f} ({accuracy*100:.2f}%)')
    print("\nClassification Report:")
    print("-" * 60)
    print(classification_rep)
    print("=" * 60)
    
    # Plot Confusion Matrix and ROC Curve
    print("\nGénération des graphiques...")
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # Confusion Matrix using ConfusionMatrixDisplay
    ConfusionMatrixDisplay.from_estimator(random_forest_model, X_test, y_test, ax=ax[0], cmap='Blues')
    ax[0].set_title('Confusion Matrix')
    
    # ROC Curve
    fpr, tpr, thresholds = roc_curve(y_test, y_pred)
    roc_auc = auc(fpr, tpr)
    
    ax[1].plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
    ax[1].plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    ax[1].set_xlim([0.0, 1.0])
    ax[1].set_ylim([0.0, 1.05])
    ax[1].set_xlabel('False Positive Rate')
    ax[1].set_ylabel('True Positive Rate')
    ax[1].set_title('Receiver Operating Characteristic (ROC) Curve')
    ax[1].legend(loc="lower right")
    
    plt.tight_layout()
    plt.savefig('ml_evaluation_plots.png', dpi=300, bbox_inches='tight')
    print("✓ Graphiques sauvegardés: ml_evaluation_plots.png")
    plt.show()


if __name__ == "__main__":
    main()
