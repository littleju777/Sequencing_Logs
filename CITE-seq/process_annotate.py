from scipy.stats import median_abs_deviation
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    
    median = np.median(M)
    MAD = median_abs_deviation(M)
    counts = outlier.sum()

    print(f"number of outliers in terms of {metric}: {counts}, median: {median}, MAD: {MAD}")
    
    return outlier


def plot_qc_feature(data):
    """
    Plot quality control features

    Parameters:
    data: annadata

    Returns:
    None
    """

    qc_info = data.obs
    
    colors = ["#cafb98", "#d1d1f6", "#ADD8E6"] 
    metrics = [
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
        ]
    ncols = 3
    nrows = 2


    fig, axes = plt.subplots(
        figsize=(15, 6), dpi=200, ncols=ncols, nrows = nrows, constrained_layout=True
    )
    
    for i in range(ncols):
        ax = axes[0,i]
        sns.violinplot(y=qc_info[metrics[i]], color=colors[i], ax=ax)
        sns.stripplot(
            y=qc_info[metrics[i]], jitter=True, color="black", ax=ax, size=0.1
        )
        ax.set_title(metrics[i])  
        
        ax = axes[1,i]
        sns.histplot(x=qc_info[metrics[i]], kde=False, color=colors[i], ax=ax)
        
    # Show the plots
    fig.tight_layout()
    plt.show()


PDAC_TME_markers = {
    'T_cells': ['Cd3d', 'Cd3e', 'Cd3g', 'Cd2', 'Cd4', 'Cd8a', 'Cd45','Foxp3'],
    'B_cells': ['Cd19', 'Cd79a', 'Cd79b', 'Ms4a1', 'Cd22', 'Ebf1', 'Ighd', 'Ighm', 'Fcmr'],
    'NK_cells': ['Ncam1', 'Klrd1', 'Nkg7', 'Gzma', 'Nkx1-1', 'Nfkbia', 'Klrb1c', 'Klrb1b', 'Klrk', 'Klra4', 'Klra7', 'Klra8', 'Klra9', 'Klre1', 'Klrg1'],
    'Monocytes': ['Cd14','Lyz2', 'Itgam', 'Ly6c2', 'Ccr2'],
    'Macrophages': ['Cd68', 'Csf1r', 'Mrc1', 'Cd86', 'Cd163', 'Itgam', 'Adgre1', 'C1qa', 'Spp1'],
    'Dendritic_cells': ['Itgax', 'Cd80', 'Itgam', 'Ccr7', 'Cd209a'],
    'Neutrophils': ['Fcgr3', 'Csf3r', 'Mpo', 'Elane', 'Ly6g', 'S100a8', 'S100a9'],
    'Epithelial_cells': ['Epcam', 'Krt18', 'Krt19', 'Cdh1', 'Muc1'],
    'Ductal_cells': ['Cftr', 'Tff1','Krt19','Sox9','Cdh1','Bmpr1a','Hnf1b','Prom1','Cldn7'],
    'Endothelial_cells': ['Cdh5','Pecam1','Cldn5'],
    'Fibrolast': ['Sfrp2','Lum','Col1a1','Pdgfra','Gja4','Cox4i2'],
    'Acinar_cells': ['Prss1', 'Klk1', 'Ctrc', 'Pnlip', 'Aldob', 'Ctrb1', 'Cela2a', 'Pnliprp1', 'Amy2b', 'Cpa1', 'Pla2g1b', 'Clps', 'Cela3b', 'Cela1', 'Prss2'],
    'Alpha_cells': ['Gcg', 'Ttr', 'Fxyd5', 'Ldb2', 'Pcsk2', 'Chga', 'Mafb', 'Pax6', 'Neurod1', 'Pyy', 'Cryba2', 'Fev', 'Irx2', 'Loxl4', 'Dpp4', 'Higd1a', 'Fap', 'Gc', 'Slc7a2', 'Gpx3', 'Tm4sf4', 'Irx2'],
    'Beta_cells': ['Nkx6-1', 'Gjd2', 'Nkx2-2', 'Slc2a2', 'Iapp', 'Ins2', 'Nkx6-1', 'Mafa', 'Ins1', 'Pdx1', 'Pax6', 'Isl1', 'Neurod1', 'Meg3', 'Pcsk2', 'G6pc2', 'Ero1lb', 'Pcsk1', 'Ffar2', 'Adcyap1', 'Dlk1', 'Pdx1', 'Hadh'],
    'Neurons': ['Map2', 'Nefl', 'Rbfox3', 'Syp', 'Tubb3']
}  

PDAC_TME_protein = {
    "B Cells": [
        "B220", "CD19", "CD20"
    ],
    "T Cells": [
        "CD2", "CD3", "CD3E", "CD4", "CD5", "CD8A", "CD8B", 
        "CD28", "CD45", "CD45.1", "CD45.2", "CD45RB", 
    ],
    "Tregs": [
        "CD25", "FOXP3", "CD127", "FR4", "CD45RA", "CD4"
    ],
    "NK Cells": [
         "CD56", "CD57", "CD160", "CD226.DNAM1", "CD244", "CD336", "NCR1", 
        "NKG2A", "NKG2C", "KLRBC_NK1.1", "NKP30", "NKP46"
    ],
    "Monocytes/Macrophages": [
        "CCR2.CD192", "CX3CR1", "F480", "CD11B", "CD14", "CD16", "CD64", "CD68", 
        "CD169", "ITAX.CD11C", "LY6C1_LY6C2", "LY6G"
    ],
    "Dendritic Cells": [
        "CD11C", "CD80", "CD86", "CD103", "CD40"
    ],
    "Neutrophils": [
        "CD16", "CD62L", "CD11B", "LY6G", "CD66b", "CXCR2", "GR1_LY6G_LY6C1_LY6C2"
    ]
}


t_markers = {
    "Activated_T" : [
        "Pdcd1", "Sell","Cd43","Cd49a","Cd49d","Klrg1","Cxcr3","Cxcr6","Cx3cr1","Cd44"
    ],
    "Exh_prog": [
        "Tcf7", "Id3", "Ccr7", "Batf", "Id3", "Sell",
        "Xcl1", "Il7r", "Slamf6", "Cxcl10", "Socs3", "Zfp36l1", "Batf", 
        "Ctla4", "Bcl6", "Xcl1", "Nr4a3", "Cxcr5", "Ccr6","Bach2","Il6ra"
    ],
    "Exh_Int": [
        "S100a6","Crip1","S100a4","Lgals3","Vim","S100a10","S100a11",
        "Lgals1","Pglyrp1"
    ],
    "Eff_like": [
         "Havcr2","Gzmb","Pdcd1","Lag3","Nr4a1", "S100a6", "Crip1", "Gzmb","Mki67",
        "S100a4", "Lgals3", "Vim", "Lsp1", "S100a10", "S100a11", "Lgals1", 
        "Pglyrp1","Birc5","Gzma","Zeb2","Bub1","Eomes"
    ],
    "Exh_klr": [
        "Klf2", "Cx3cr1", "Zeb2", "Klre1", "Klrc1", "Klrd1", 
        "Klrg1", "Klra3"
    ],
    "Exh_term" : [
    "Cd7", "Cd160", "Rgs1", "Cxcr6",  
     "Nr4a2", "Tigit", "Lag3", "Sh2d2a", 
    "Glrx",  "Pdcd1",  "Trac",  "Hic1", "Id2", "Bhlhe41", 
     "Cxcr3", "Plac8", "Tox", "Prdm1"
    ]
 
}


tam_markers = {
    "C1qa+ Mrc1+ Adgre1+ TAM" : [
    'Ccl8', 'Ms4a7', 'C1qa', 'C1qc', 'Cd63',
    'C1qb', 'Cd81', 'Sepp1', 'Ccl12', 'Apoe',
    'Vcam1', 'Pf4', 'Cd72', 'Timp2', 'Mrc1'
    ],
    "Arg1+ Vegfa+ TAM" : [
    'Arg1', 'Cxcl1', 'Spp1', 'Mmp12', 'Hilpda',
    'Cxcl2', 'Ero1l', 'Cxcl3', 'Hmox1', 'Bnip3',
    'Ccl2', 'Ctsl', 'Pf4', 'Ndrg1'
    ],
    "Profliferating TAM" : [
    "Top2a", "Mki67"
    ],
    "ISG Mono/Macro" : [
    "Ifit3", "Isg15"
    ],
    "AP Monocyte" : [
    "H2-DMb1","H2-DMa"
    ],
    "Chil3+ Hp+ Monocyte" : [
    'Chil3', 'Hp', 'Cxcl10', 'Ifitm6', 'Gsr',
    'Plac8', 'Mgst1', 'Samhd1', 'Slpi', 'Smpdl3a',
    'Ifi47', 'Sell'
    ],
    "Non-classical Monocyte" : [
    'Ear2', 'Pglyrp1', 'Fabp4', 'Ace', 'Gngt2',
    'Ccl3', 'Treml4', 'Nr4a1', 'Eno3', 'Itgal',
    'Adgre4', 'Spic', 'Stap1', 'Cd36', 'Gsr'
    ],
    "TAM" : [
    "Ly6c2", "Ccr2", "Adgre1"
    ]



}



