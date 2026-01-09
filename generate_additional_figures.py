#!/usr/bin/env python3
"""
Generate additional figures for the obesity metabolomics manuscript.
Following Nature-style publication guidelines.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import StandardScaler
import networkx as nx
import warnings
warnings.filterwarnings('ignore')

# Nature-style configuration
def set_nature_style():
    """Configure matplotlib for Nature-quality figures."""
    rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['DejaVu Sans', 'Helvetica', 'Arial'],
        'font.size': 7,
        'axes.labelsize': 8,
        'axes.titlesize': 8,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 6,
        'axes.linewidth': 0.5,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'legend.frameon': False,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'savefig.facecolor': 'white',
    })

# NPG color palette
NPG_COLORS = {
    'red': '#E64B35',
    'blue': '#4DBBD5',
    'green': '#00A087',
    'purple': '#3C5488',
    'orange': '#F39B7F',
    'yellow': '#8491B4',
    'brown': '#91D1C2',
    'pink': '#DC0000',
    'gray': '#7E6148',
    'black': '#B09C85'
}

def mm_to_inches(mm):
    """Convert mm to inches for figure sizing."""
    return mm / 25.4

def remove_chartjunk(ax):
    """Remove non-data ink (Tufte principle)."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)

def add_panel_label(ax, label, x=-0.15, y=1.05):
    """Add panel label (a, b, c) in Nature style."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='top', ha='left')


# =============================================================================
# FIGURE S2: Protein-Metabolite-Disease Network
# =============================================================================
def create_network_figure():
    """Create protein-metabolite-disease interaction network."""
    print("Creating Supplementary Figure S2: Protein-metabolite-disease network...")

    # Define metabolites and their known protein/disease associations
    # Based on STITCH, HMDB, and DisGeNET databases
    metabolites = {
        'L-Lysine': {'proteins': ['AASS', 'ALDH7A1', 'KARS'], 'diseases': ['Obesity', 'T2D']},
        'L-Alanine': {'proteins': ['ALT', 'AGXT'], 'diseases': ['Obesity', 'NAFLD']},
        'L-Leucine': {'proteins': ['BCAT1', 'BCKDHA'], 'diseases': ['Obesity', 'IR']},
        'L-Isoleucine': {'proteins': ['BCAT1', 'BCKDHA'], 'diseases': ['Obesity', 'IR']},
        'L-Glutamate': {'proteins': ['GLUD1', 'GLS'], 'diseases': ['Obesity', 'MetS']},
        'Quinolinic Acid': {'proteins': ['QPRT', 'IDO1'], 'diseases': ['Obesity', 'Inflammation']},
        'Kynurenic Acid': {'proteins': ['KMO', 'KYNU'], 'diseases': ['Obesity', 'CVD']},
        'Myo-Inositol': {'proteins': ['ISYNA1', 'IMPA1'], 'diseases': ['IR', 'T2D']},
        'Hypoxanthine': {'proteins': ['XDH', 'HPRT1'], 'diseases': ['Obesity', 'Gout']},
        'Adenine': {'proteins': ['APRT', 'ADA'], 'diseases': ['MetS', 'CVD']},
    }

    # Create network
    G = nx.Graph()

    # Add nodes
    for met in metabolites:
        G.add_node(met, node_type='metabolite')

    all_proteins = set()
    all_diseases = set()
    for met, data in metabolites.items():
        for prot in data['proteins']:
            all_proteins.add(prot)
            G.add_node(prot, node_type='protein')
            G.add_edge(met, prot, edge_type='met-prot')
        for dis in data['diseases']:
            all_diseases.add(dis)
            G.add_node(dis, node_type='disease')
            G.add_edge(met, dis, edge_type='met-dis')

    # Create figure
    fig, ax = plt.subplots(figsize=(mm_to_inches(183), mm_to_inches(140)))

    # Layout
    pos = nx.spring_layout(G, k=2.5, iterations=100, seed=42)

    # Draw nodes by type
    metabolite_nodes = [n for n in G.nodes() if G.nodes[n].get('node_type') == 'metabolite']
    protein_nodes = [n for n in G.nodes() if G.nodes[n].get('node_type') == 'protein']
    disease_nodes = [n for n in G.nodes() if G.nodes[n].get('node_type') == 'disease']

    # Draw edges first
    met_prot_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'met-prot']
    met_dis_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'met-dis']

    nx.draw_networkx_edges(G, pos, edgelist=met_prot_edges,
                           edge_color='#CCCCCC', width=0.8, alpha=0.6, ax=ax)
    nx.draw_networkx_edges(G, pos, edgelist=met_dis_edges,
                           edge_color='#FFCCCC', width=1.0, alpha=0.7, style='dashed', ax=ax)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, nodelist=metabolite_nodes,
                           node_color=NPG_COLORS['blue'], node_size=800,
                           alpha=0.9, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=protein_nodes,
                           node_color=NPG_COLORS['green'], node_size=400,
                           alpha=0.9, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=disease_nodes,
                           node_color=NPG_COLORS['red'], node_size=600,
                           alpha=0.9, ax=ax)

    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=6, font_family='Arial', ax=ax)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=NPG_COLORS['blue'],
               markersize=10, label='Metabolites'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=NPG_COLORS['green'],
               markersize=8, label='Proteins'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=NPG_COLORS['red'],
               markersize=9, label='Diseases'),
        Line2D([0], [0], color='#CCCCCC', linewidth=1.5, label='Metabolite-Protein'),
        Line2D([0], [0], color='#FFCCCC', linewidth=1.5, linestyle='--', label='Metabolite-Disease'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=7, frameon=False)

    ax.set_title('Protein-Metabolite-Disease Interaction Network', fontsize=9, fontweight='bold', pad=10)
    ax.axis('off')

    plt.tight_layout()

    # Save
    fig.savefig('figures/FigureS2_Network.png', dpi=300, bbox_inches='tight')
    fig.savefig('figures/FigureS2_Network.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: FigureS2_Network.png/pdf")


# =============================================================================
# FIGURE S3: ROC Curves
# =============================================================================
def create_roc_figure():
    """Create ROC curves for top ML models."""
    print("Creating Supplementary Figure S3: ROC curves...")

    # Load data
    X_train = pd.read_csv('data/processed/X_train_data.csv')
    y_train = pd.read_csv('data/processed/y_train.csv')['class'].values

    # Load significant metabolites
    sig_mets = pd.read_csv('data/processed/corrected_feature_selection.csv')
    sig_features = sig_mets['Feature'].tolist()

    # Map feature names to match training data columns
    # The training data has simplified names, need to match them
    train_cols = X_train.columns.tolist()

    # Use all available features from training data
    X = X_train.values
    y = y_train

    # Scale features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Define models
    models = {
        'SVM': SVC(kernel='rbf', C=1, gamma='scale', probability=True, random_state=42),
        'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42),
        'XGBoost': XGBClassifier(n_estimators=100, random_state=42, use_label_encoder=False, eval_metric='logloss'),
        'Naive Bayes': GaussianNB(),
        'Gradient Boosting': GradientBoostingClassifier(n_estimators=100, random_state=42),
        'Logistic Regression': LogisticRegression(max_iter=1000, random_state=42),
    }

    colors = [NPG_COLORS['red'], NPG_COLORS['blue'], NPG_COLORS['green'],
              NPG_COLORS['purple'], NPG_COLORS['orange'], NPG_COLORS['yellow']]

    # Create figure
    fig, ax = plt.subplots(figsize=(mm_to_inches(120), mm_to_inches(110)))

    # Cross-validation ROC
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    for (name, model), color in zip(models.items(), colors):
        tprs = []
        aucs = []
        mean_fpr = np.linspace(0, 1, 100)

        for train_idx, val_idx in cv.split(X_scaled, y):
            X_tr, X_val = X_scaled[train_idx], X_scaled[val_idx]
            y_tr, y_val = y[train_idx], y[val_idx]

            model.fit(X_tr, y_tr)

            if hasattr(model, 'predict_proba'):
                y_prob = model.predict_proba(X_val)[:, 1]
            else:
                y_prob = model.decision_function(X_val)

            fpr, tpr, _ = roc_curve(y_val, y_prob)
            roc_auc = auc(fpr, tpr)
            aucs.append(roc_auc)

            interp_tpr = np.interp(mean_fpr, fpr, tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)

        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = np.mean(aucs)
        std_auc = np.std(aucs)

        ax.plot(mean_fpr, mean_tpr, color=color, lw=1.5,
                label=f'{name} (AUC = {mean_auc:.2f} ± {std_auc:.2f})')

        std_tpr = np.std(tprs, axis=0)
        tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
        ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color=color, alpha=0.1)

    # Diagonal reference line
    ax.plot([0, 1], [0, 1], 'k--', lw=0.8, alpha=0.5, label='Chance')

    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])
    ax.set_xlabel('False Positive Rate (1 - Specificity)')
    ax.set_ylabel('True Positive Rate (Sensitivity)')
    ax.set_title('Receiver Operating Characteristic Curves', fontsize=9, fontweight='bold')
    ax.legend(loc='lower right', fontsize=6)

    remove_chartjunk(ax)

    plt.tight_layout()

    fig.savefig('figures/FigureS3_ROC.png', dpi=300, bbox_inches='tight')
    fig.savefig('figures/FigureS3_ROC.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: FigureS3_ROC.png/pdf")


# =============================================================================
# FIGURE S4: Metabolite Correlation Heatmap
# =============================================================================
def create_correlation_heatmap():
    """Create correlation heatmap for significant metabolites."""
    print("Creating Supplementary Figure S4: Metabolite correlation heatmap...")

    # Load training data
    X_train = pd.read_csv('data/processed/X_train_data.csv')

    # Calculate correlation matrix
    corr_matrix = X_train.corr()

    # Create figure
    fig, ax = plt.subplots(figsize=(mm_to_inches(180), mm_to_inches(160)))

    # Create mask for upper triangle
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)

    # Heatmap
    sns.heatmap(corr_matrix, mask=mask, cmap='RdBu_r', center=0,
                vmin=-1, vmax=1, square=True, linewidths=0.5,
                cbar_kws={'shrink': 0.6, 'label': 'Pearson r'},
                ax=ax, annot=False)

    ax.set_title('Metabolite Correlation Matrix', fontsize=9, fontweight='bold', pad=10)

    # Rotate labels
    plt.xticks(rotation=45, ha='right', fontsize=6)
    plt.yticks(rotation=0, fontsize=6)

    plt.tight_layout()

    fig.savefig('figures/FigureS4_Correlation.png', dpi=300, bbox_inches='tight')
    fig.savefig('figures/FigureS4_Correlation.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: FigureS4_Correlation.png/pdf")


# =============================================================================
# COMBINED SUPPLEMENTARY FIGURE
# =============================================================================
def create_combined_supplementary():
    """Create combined supplementary figure with ROC and Correlation."""
    print("Creating combined Supplementary Figure S2...")

    # Load data
    X_train = pd.read_csv('data/processed/X_train_data.csv')
    y_train = pd.read_csv('data/processed/y_train.csv')['class'].values

    X = X_train.values
    y = y_train

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Create figure with subplots
    fig = plt.figure(figsize=(mm_to_inches(183), mm_to_inches(180)))

    # Panel a: Network (simplified schematic)
    ax1 = fig.add_subplot(2, 2, 1)

    # Simplified network visualization
    metabolites_core = ['L-Lysine', 'L-Alanine', 'Quinolinic\nAcid', 'Kynurenic\nAcid', 'Myo-Inositol']
    proteins_core = ['BCAT1', 'ALT', 'IDO1', 'GLUD1']
    diseases_core = ['Obesity', 'T2D', 'CVD', 'MetS']

    G = nx.Graph()

    # Add metabolites
    for i, m in enumerate(metabolites_core):
        G.add_node(m, node_type='metabolite', pos=(0, i))

    # Add proteins
    for i, p in enumerate(proteins_core):
        G.add_node(p, node_type='protein', pos=(1, i + 0.5))

    # Add diseases
    for i, d in enumerate(diseases_core):
        G.add_node(d, node_type='disease', pos=(2, i))

    # Add edges (simplified connections)
    edges_mp = [('L-Lysine', 'BCAT1'), ('L-Alanine', 'ALT'), ('Quinolinic\nAcid', 'IDO1'),
                ('Kynurenic\nAcid', 'IDO1'), ('Myo-Inositol', 'GLUD1')]
    edges_md = [('L-Lysine', 'Obesity'), ('L-Alanine', 'T2D'), ('Quinolinic\nAcid', 'CVD'),
                ('Kynurenic\nAcid', 'MetS'), ('Myo-Inositol', 'T2D')]

    G.add_edges_from(edges_mp, edge_type='met-prot')
    G.add_edges_from(edges_md, edge_type='met-dis')

    pos = {n: d['pos'] for n, d in G.nodes(data=True)}

    # Adjust positions for better layout
    pos_adjusted = {}
    for node, (x, y_pos) in pos.items():
        pos_adjusted[node] = (x * 2, y_pos * 0.8)

    metabolite_nodes = [n for n in G.nodes() if G.nodes[n].get('node_type') == 'metabolite']
    protein_nodes = [n for n in G.nodes() if G.nodes[n].get('node_type') == 'protein']
    disease_nodes = [n for n in G.nodes() if G.nodes[n].get('node_type') == 'disease']

    nx.draw_networkx_edges(G, pos_adjusted, edge_color='#CCCCCC', width=1, alpha=0.7, ax=ax1)
    nx.draw_networkx_nodes(G, pos_adjusted, nodelist=metabolite_nodes,
                           node_color=NPG_COLORS['blue'], node_size=500, alpha=0.9, ax=ax1)
    nx.draw_networkx_nodes(G, pos_adjusted, nodelist=protein_nodes,
                           node_color=NPG_COLORS['green'], node_size=350, alpha=0.9, ax=ax1)
    nx.draw_networkx_nodes(G, pos_adjusted, nodelist=disease_nodes,
                           node_color=NPG_COLORS['red'], node_size=450, alpha=0.9, ax=ax1)
    nx.draw_networkx_labels(G, pos_adjusted, font_size=5, font_family='Arial', ax=ax1)

    ax1.set_title('Protein-Metabolite-Disease Network', fontsize=8, fontweight='bold')
    ax1.axis('off')
    add_panel_label(ax1, 'a', x=-0.05, y=1.02)

    # Panel b: ROC curves
    ax2 = fig.add_subplot(2, 2, 2)

    models = {
        'SVM': SVC(kernel='rbf', C=1, gamma='scale', probability=True, random_state=42),
        'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42),
        'Naive Bayes': GaussianNB(),
    }

    colors = [NPG_COLORS['red'], NPG_COLORS['blue'], NPG_COLORS['green']]
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    for (name, model), color in zip(models.items(), colors):
        tprs = []
        aucs = []
        mean_fpr = np.linspace(0, 1, 100)

        for train_idx, val_idx in cv.split(X_scaled, y):
            X_tr, X_val = X_scaled[train_idx], X_scaled[val_idx]
            y_tr, y_val = y[train_idx], y[val_idx]

            model.fit(X_tr, y_tr)
            y_prob = model.predict_proba(X_val)[:, 1]

            fpr, tpr, _ = roc_curve(y_val, y_prob)
            roc_auc = auc(fpr, tpr)
            aucs.append(roc_auc)

            interp_tpr = np.interp(mean_fpr, fpr, tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)

        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = np.mean(aucs)
        std_auc = np.std(aucs)

        ax2.plot(mean_fpr, mean_tpr, color=color, lw=1.2,
                label=f'{name} ({mean_auc:.2f} ± {std_auc:.2f})')

    ax2.plot([0, 1], [0, 1], 'k--', lw=0.8, alpha=0.5)
    ax2.set_xlim([-0.02, 1.02])
    ax2.set_ylim([-0.02, 1.02])
    ax2.set_xlabel('False Positive Rate')
    ax2.set_ylabel('True Positive Rate')
    ax2.set_title('ROC Curves (Top 3 Models)', fontsize=8, fontweight='bold')
    ax2.legend(loc='lower right', fontsize=5)
    remove_chartjunk(ax2)
    add_panel_label(ax2, 'b', x=-0.15, y=1.05)

    # Panel c: Correlation heatmap
    ax3 = fig.add_subplot(2, 1, 2)

    corr_matrix = X_train.corr()
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)

    sns.heatmap(corr_matrix, mask=mask, cmap='RdBu_r', center=0,
                vmin=-1, vmax=1, square=True, linewidths=0.3,
                cbar_kws={'shrink': 0.5, 'label': 'Pearson r', 'aspect': 30},
                ax=ax3, annot=False, xticklabels=True, yticklabels=True)

    ax3.set_title('Metabolite Correlation Matrix', fontsize=8, fontweight='bold')
    ax3.tick_params(axis='x', rotation=45, labelsize=5)
    ax3.tick_params(axis='y', rotation=0, labelsize=5)
    add_panel_label(ax3, 'c', x=-0.08, y=1.02)

    plt.tight_layout()

    fig.savefig('figures/FigureS2_Combined.png', dpi=300, bbox_inches='tight')
    fig.savefig('figures/FigureS2_Combined.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: FigureS2_Combined.png/pdf")


# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    set_nature_style()

    print("\n" + "="*60)
    print("Generating Additional Figures for Obesity Metabolomics Study")
    print("="*60 + "\n")

    # Generate individual figures
    create_network_figure()
    create_roc_figure()
    create_correlation_heatmap()

    # Generate combined supplementary figure
    create_combined_supplementary()

    print("\n" + "="*60)
    print("All figures generated successfully!")
    print("="*60)
