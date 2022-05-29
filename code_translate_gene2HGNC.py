# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 19:05:19 2021

@author: da.salazarb
"""

import cobra
import pandas as pd

## genes a mapear
# GSE10142 = pd.read_excel("D:/SB-HCC/Regulons_intersects .xlsx", sheet_name="GSE10142")
GSE14520 = pd.read_excel("D:/SB-HCC/Regulons_intersects .xlsx", sheet_name="GSE14520")

GSE14520_AEBP1 = GSE14520.AEBP1.loc[~GSE14520.AEBP1.isna().values,].str.split(pat="///").explode("AEBP1")
GSE14520_AR = GSE14520.AR.loc[~GSE14520.AR.isna().values,].str.split(pat="///").explode("AR")
GSE14520_NR1I3 = GSE14520.NR1I3.loc[~GSE14520.NR1I3.isna().values,].str.split(pat="///").explode("NR1I3")

salida = pd.concat([GSE14520_AEBP1, GSE14520_AR, GSE14520_NR1I3], axis=1)
salida.to_csv("D:/SB-HCC/newGSE14520.csv")

#### 
GSE10142_AEBP1 = pd.read_table("D:/SB-HCC/GSE10142_AEBP1.csv",sep=";"); GSE10142_AEBP1_v2 = list(GSE10142_AEBP1.loc[:,"HGNC ID"].values)
GSE10142_AR = pd.read_table("D:/SB-HCC/GSE10142_AR.csv", sep=";"); GSE10142_AR_v2 = list(GSE10142_AR.loc[:,"HGNC ID"].values)
GSE10142_NR1I3 = pd.read_table("D:/SB-HCC/GSE10142_NR1I3.csv", sep=";"); GSE10142_NR1I3_v2 = list(GSE10142_NR1I3.loc[:,"HGNC ID"].values)

#### 
GSE14520_AEBP1 = pd.read_table("D:/SB-HCC/GSE14520_AEBP1.csv", sep=";"); GSE14520_AEBP1_v2 = list(GSE14520_AEBP1.loc[:,"HGNC ID"].values)
GSE14520_AR = pd.read_table("D:/SB-HCC/GSE14520_AR.csv", sep=";"); GSE14520_AR_v2 = list(GSE14520_AR.loc[:,"HGNC ID"].values)
GSE14520_NR1I3 = pd.read_table("D:/SB-HCC/GSE14520_NR1I3.csv", sep=";"); GSE14520_NR1I3_v2 = list(GSE14520_NR1I3.loc[:,"HGNC ID"].values)

## modelo
model = cobra.io.read_sbml_model('D:/SB-HCC/MODEL1603150001.xml')

# %% 
genes_lista_recon = []
for rxn in model.reactions:
    rxn_temp = list(rxn.genes)
    for gene in rxn_temp:
        genes_lista_recon.append(gene.id)
        
genes_lista_recon = list(set(genes_lista_recon))

# %%
GSE10142_AEBP1_v2 = list(set(genes_lista_recon) & set(GSE10142_AEBP1_v2))
GSE10142_AR_v2 = list(set(genes_lista_recon) & set(GSE10142_AR_v2))
GSE10142_NR1I3_v2 = list(set(genes_lista_recon) & set(GSE10142_NR1I3_v2))

GSE14520_AEBP1_v2 = list(set(genes_lista_recon) & set(GSE14520_AEBP1_v2))
GSE14520_AR_v2 = list(set(genes_lista_recon) & set(GSE14520_AR_v2))
GSE14520_NR1I3_v2 = list(set(genes_lista_recon) & set(GSE14520_NR1I3_v2))

### hgnc to genesymbol
GSE10142_AEBP1 = GSE10142_AEBP1.Input[GSE10142_AEBP1["HGNC ID"].isin(GSE10142_AEBP1_v2)].reset_index(drop = True)
GSE10142_AR = GSE10142_AR.Input[GSE10142_AR["HGNC ID"].isin(GSE10142_AR_v2)].reset_index(drop = True)
GSE10142_NR1I3 = GSE10142_NR1I3.Input[GSE10142_NR1I3["HGNC ID"].isin(GSE10142_NR1I3_v2)].reset_index(drop = True)

GSE14520_AEBP1 = GSE14520_AEBP1.Input[GSE14520_AEBP1["HGNC ID"].isin(GSE14520_AEBP1_v2)].reset_index(drop = True)
GSE14520_AR = GSE14520_AR.Input[GSE14520_AR["HGNC ID"].isin(GSE14520_AR_v2)].reset_index(drop = True)
GSE14520_NR1I3 = GSE14520_NR1I3.Input[GSE14520_NR1I3["HGNC ID"].isin(GSE14520_NR1I3_v2)].reset_index(drop = True)

salida = pd.concat([GSE10142_AEBP1,GSE10142_AR, GSE10142_NR1I3], axis=1)
salida.to_csv("D:/SB-HCC/mappedGSE10142.csv")

salida = pd.concat([GSE14520_AEBP1,GSE14520_AR, GSE14520_NR1I3], axis=1)
salida.to_csv("D:/SB-HCC/mappedGSE14520.csv")
