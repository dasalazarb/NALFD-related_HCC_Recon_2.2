# -*- coding: utf-8 -*-
"""
Created on Sun Oct 31 15:08:53 2021

@author: da.salazarb
"""
# %% Libraries and data
import cobra
import pandas as pd
import numpy as np

GSE14520_MoA = pd.read_table("D:/SB-HCC/00_MoA/GSE14520_MoA.txt")
GSE10142_MoA = pd.read_table("D:/SB-HCC/00_MoA/GSE10142_MoA.txt")

neg_AEBP1 = list(GSE14520_MoA.loc[(GSE14520_MoA["modeOfAction"] == -1) & (GSE14520_MoA["from"] == "AEBP1"),:]["to_hgnc"])
neg_AR = list(GSE14520_MoA.loc[(GSE14520_MoA["modeOfAction"] == -1) & (GSE14520_MoA["from"] == "AR"),:]["to_hgnc"])
neg_NR1I3 = list(GSE14520_MoA.loc[(GSE14520_MoA["modeOfAction"] == -1) & (GSE14520_MoA["from"] == "NR1I3"),:]["to_hgnc"])

pos_AEBP1 = list(GSE14520_MoA.loc[(GSE14520_MoA["modeOfAction"] == 1) & (GSE14520_MoA["from"] == "AEBP1"),:]["to_hgnc"])
pos_AR = list(GSE14520_MoA.loc[(GSE14520_MoA["modeOfAction"] == 1) & (GSE14520_MoA["from"] == "AR"),:]["to_hgnc"])
pos_NR1I3 = list(GSE14520_MoA.loc[(GSE14520_MoA["modeOfAction"] == 1) & (GSE14520_MoA["from"] == "NR1I3"),:]["to_hgnc"])

## modelo
model = cobra.io.read_sbml_model('D:/SB-HCC/MODEL1603150001.xml')

pos_neg = {"neg_AEBP1":neg_AEBP1,"neg_AR":neg_AR, "neg_NR1I3":neg_NR1I3, 
           "pos_AEBP1":pos_AEBP1, "pos_AR":pos_AR, "pos_NR1I3": pos_NR1I3}

# %% fucntions
def get_list_rxn(rxn_temp,pos_temp,genes,genes_rxn):
    if len(set(rxn_temp).intersection(pos_temp)) > 0: ## hay genes para revisar?
        indices = []
        list_rxn_names_temp = []
        for j in set(rxn_temp).intersection(pos_temp): ## recorre genes para apagar o prender
            indices = indices + list(np.where(j == genes)[0])
        
        indices = np.array(indices)
        
        for k in genes_rxn[indices]: ## para cada gen revisar las reacciones
            list_rxn_names_temp = [i.name for i in list(k.reactions)]
    else:
        list_rxn_names_temp = []
        
    return list_rxn_names_temp
    
def modify_model(model, pos_neg):
    
    genes = np.array([i.name for i in model.genes]) # nombres genes
    rxns  = np.array([i.name for i in model.reactions]) # nombres reacciones
    genes_rxn = pd.DataFrame(model.genes).iloc[:,0] # identificador gen
    
    pos_list = []
    neg_list = []
    for rxn in model.reactions: ## recorre reacciones
        rxn_temp = [gene.id for gene in rxn.genes]
        
        ### positive MoA
        pos_temp = {k: v for k,v in pos_neg.items() if "pos" in k}
        pos_temp = [i for v in pos_temp.values() for i in v]
        
        pos_list = pos_list + get_list_rxn(rxn_temp,pos_temp,genes,genes_rxn)
        
        ### negative MoA
        neg_temp = {k: v for k,v in pos_neg.items() if "neg" in k}
        neg_temp = [i for v in neg_temp.values() for i in v]
        
        neg_list = neg_list + get_list_rxn(rxn_temp,neg_temp,genes,genes_rxn)
        
    pos_list = list(np.unique(pos_list))
    neg_list = list(np.unique(neg_list))
    
    lista_indices_pos = []
    for l in pos_list: ## modificar lb y ub para cada reaccion
        indice_rxn = list(np.where(l == rxns)[0])
        lista_indices_pos = lista_indices_pos + indice_rxn
        for indice in indice_rxn:
            model.reactions[indice].upper_bound = 1000
            model.reactions[indice].lower_bound = -1000
        
    lista_indices_neg = []
    for l in neg_list: ## modificar lb y ub para cada reaccion
        indice_rxn = list(np.where(l == rxns)[0])
        lista_indices_neg = lista_indices_neg + indice_rxn
        for indice in indice_rxn:
            # print(indice)
            # print(model.reactions[indice].name)
            model.reactions[indice].lower_bound = 0
            model.reactions[indice].upper_bound = 10
            # print(model.reactions[indice].lower_bound)
            # print(model.reactions[indice].upper_bound)
            
    return model, lista_indices_pos, lista_indices_neg

# %% to generate the model of NAFLD-related HCC model
model_temp, lista_indices_pos, lista_indices_neg = modify_model(model, pos_neg)

## Check point
## positive genes
rxn_temp = model_temp.reactions[lista_indices_pos[10]]
rxn_temp.name
rxn_temp.reaction
rxn_temp.lower_bound
rxn_temp.upper_bound

dict_pos = dict()
dict_pos["reaction"] = []; dict_pos["gene"] = []; dict_pos["name"] = []
for i in lista_indices_pos:
    rxn_temp = model_temp.reactions[i]
    dict_pos["name"].append(rxn_temp.name)
    dict_pos["reaction"].append(rxn_temp.reaction)
    dict_pos["gene"].append(rxn_temp.gene_name_reaction_rule)

## negative gene
rxn_temp = model_temp.reactions[lista_indices_neg[0]]
rxn_temp.name
rxn_temp.reaction
rxn_temp.lower_bound
rxn_temp.upper_bound

dict_neg = dict()
dict_neg["reaction"] = []; dict_neg["gene"] = []; dict_neg["name"] = []
for i in lista_indices_neg:
    rxn_temp = model_temp.reactions[i]
    dict_neg["name"].append(rxn_temp.name)
    dict_neg["reaction"].append(rxn_temp.reaction)
    dict_neg["gene"].append(rxn_temp.gene_name_reaction_rule)

## which reactions where modified in NALFD-related HCC
df_pos = pd.DataFrame.from_dict(dict_pos); df_pos.to_csv("D:/SB-HCC/which_genes_regulated_pos.csv", sep=";")
df_neg = pd.DataFrame.from_dict(dict_neg); df_neg.to_csv("D:/SB-HCC/which_genes_regulated_neg.csv", sep=";")

# %% run pFBA
model.objective = 'ATPS4m'
# fba_solution = model.optimize()
pfba_solution = cobra.flux_analysis.pfba(model)

model_temp.objective = 'ATPS4m'
# fba_solution = model_temp.optimize()
pfba_solution = cobra.flux_analysis.pfba(model_temp)
pfba_solution.status
