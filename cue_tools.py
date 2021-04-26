import tempfile
import cobra
import requests
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def loadModelURL(url, fmt):
    """
    Load a cobra model from a URL
    Inputs:
    | url <str>: URL to metabolic model
    | fmt <str>: Format of metabolic model 
        options: {'sbml', 'json', 'mat'}
    Outputs:
    | model: cobra model
    
    Example: 
    1)
    url = 'http://bigg.ucsd.edu/static/models/e_coli_core.mat'
    model = loadModelURL(url, 'mat')
    
    2)
    url = 'http://bigg.ucsd.edu/static/models/e_coli_core.xml'
    model = loadModelURL(url, 'sbml')
    """
    # Select load function
    load_functions = {'sbml': cobra.io.read_sbml_model,
                      'json': cobra.io.load_json_model,
                      'mat': cobra.io.load_matlab_model}
    if fmt not in load_functions:
        raise KeyError(f'{fmt} is not a valid format. Choose from {set(load_functions.keys())}')
    load_function = load_functions[fmt]
    
    # Load model from URL
    with tempfile.NamedTemporaryFile(delete=False) as tf:
        resp = requests.get(url)
        tf.write(resp.content)
        model = load_function(tf.name)
    os.remove(tf.name)
    return model

def computeCUE(model, biomass_rxn, co2_rxn):
    """
    Compute carbon use efficiency = rate of biomass production / rate of CO2 secretion
    Inputs:
    | model: Cobra model
    | {biomass, co2}_rxn: Model ID of biomass and co2 secretion reaction
    Outputs:
    | cue: biomass rate / co2 rate
    
    Example:
    # Load model
    url = 'http://bigg.ucsd.edu/static/models/e_coli_core.mat'
    model = loadModelURL(url, 'mat')
    # Calculate cue
    cue = computeCUE(model, 'BIOMASS_Ecoli_core_w_GAM', 'EX_co2_e')
    """
    # Run FBA
    res = model.optimize()
    # Compute CUE
    cue = res.get_primal_by_id(biomass_rxn)/res.get_primal_by_id(co2_rxn)
    return cue

def computeEpistasis (model,objective_func = "default",heatmap = False, labels = False, export_matrix = False):
    """
    Function for computing epistatic reactions for a metabolic model
    Does not ouput values of single knockouts
    Inputs:
    | model: Cobra model
    | objective_func <str>: Objective function, default is model default
    | heatmap <bool>: True of false, true generates a heatmap, default is false
    | labels <bool>: True or false, true add labels to the heatmap
    | export_matrix <bool>: True of false, true exports heatmap matrix
    Outputs:
    | epi_dist: Distribution of epistatic interactions 
    | heatmap of epistatic interactions
    
    Example:
    # To compute distribution of epistatic interactions
    epi = computeEpistasis(model)
    # To compute the distribution and plot the heatmap
    epi = computeEpistasis(model, heatmap = True)
    """
    #set model objective function
    if objective_func != "default":
        model.objective = [objective_func]
        
    #compute wild type growth rate
    solution = model.optimize()
    wt_grow = solution.objective_value
    
    ## knock-outs loops ##
    # initialize empty matricies
    rxns = len(model.reactions)
    single_ko = np.zeros((rxns))
    v1v2_grow = np.zeros((rxns,rxns))

    ## Single knockouts ##
    for i in range(rxns):
        #buffer
        upper_i = model.reactions[i].upper_bound
        lower_i = model.reactions[i].lower_bound
        # set upper and lower bounds to zero
        model.reactions[i].upper_bound = 0
        model.reactions[i].lower_bound = 0
        # solve model and record growth rate
        solution = model.optimize()
        single_ko[i] = solution.objective_value
        # return bounds to their previous state
        model.reactions[i].upper_bound = upper_i
        model.reactions[i].lower_bound = lower_i
    
    ## Combo knockout ##
    for i in range(rxns):
        for j in range(rxns):
            if j > i:
                #buffer
                upper_i = model.reactions[i].upper_bound
                lower_i = model.reactions[i].lower_bound
                upper_j = model.reactions[j].upper_bound
                lower_j = model.reactions[j].lower_bound
                # Set bounds on rxns to zero, now both are zero
                model.reactions[i].upper_bound = 0
                model.reactions[i].lower_bound = 0
                model.reactions[j].upper_bound = 0
                model.reactions[j].lower_bound = 0
                # Solve
                solution = model.optimize()
                v1v2_grow[i,j] = solution.objective_value
                # return them to what they were
                model.reactions[i].upper_bound = upper_i
                model.reactions[i].lower_bound = lower_i
                model.reactions[j].upper_bound = upper_j
                model.reactions[j].lower_bound = lower_j

    #adjusting matricies 
    v1_grow = single_ko + np.zeros((rxns,rxns))
    v2_grow = np.transpose(v1_grow)
    np.fill_diagonal(v1v2_grow,single_ko)
    
    # distribution of epistatic interactions
    epistasis = (v1v2_grow/wt_grow) - ((v1_grow/wt_grow) * (v2_grow/wt_grow))
    ep_dist = epistasis[np.triu_indices(rxns,1)]
    epistasis_full = np.triu(epistasis)+np.rot90(np.fliplr(np.triu(epistasis,1)))
            
    # heatmap
    if heatmap:
        plt.figure(figsize = (12,9))
        if labels:
            reactions = []
            for x in model.reactions:
                reactions.append(x.id)
            sns.heatmap(epistasis_full,xticklabels = reactions, yticklabels = reactions, cmap='mako', linecolor = 'dimgrey', linewidth = 0.005)
        else:
            sns.heatmap(epistasis_full, cmap='mako', linecolor = 'dimgrey', linewidth = 0.005)
        plt.show
    
    #export distribution and matrix
    if export_matrix:
        return(ep_dist,epistasis_full)
    
    return(ep_dist)

def CUEDef(model):
    
    
    # Getting the Uptake and Secretion fluxes from the model 
    Summary = model.summary()
    UptakeSummary = Summary._display_flux(Summary.uptake_flux, False, 'C', 1e-7)
    SecretionSummary = Summary._display_flux(Summary.secretion_flux, False, 'C', 1e-7)
    
    
    # First CUE Definition
    
    # Finding What Uptake Fluxes Have Metabolites with Carbon

    UFC = UptakeSummary[UptakeSummary['C-Number'] > 0].index.values

    # Finding What Secretion Fluxes Have Metabolites with Carbon

    SFC = SecretionSummary[SecretionSummary['C-Number']>0].index.values

    CUE1 = 0
    CUE2 = 0
    for n in UFC:

        CUE1 = CUE1 + UptakeSummary.loc[n,'flux']*UptakeSummary.loc[n,'C-Number']

    for n in SFC:

        CUE2 = CUE2 + SecretionSummary.loc[n,'flux']*SecretionSummary.loc[n,'C-Number']
        
    FinalCUE1 = (CUE1 + CUE2)/CUE1
    
    # Second CUE Definition
    
    for n in SecretionSummary.metabolite.values:
    
        Forms = getattr(model.metabolites,n).formula
    
        if (Forms == 'CO2'):
        
            # Find the location of the flux for Co2
        
            L = SecretionSummary[SecretionSummary['metabolite']==n].index.values
        
            CUE22 = SecretionSummary.loc[L,'flux'].values
        
    FinalCUE2 = (CUE1 + CUE22)/CUE1
    
    # Third CUE Definition
    
    solution = model.optimize()
    BiomassProd = solution.objective_value
    CUE23 = -CUE2 + BiomassProd
    FinalCUE3 = BiomassProd/(CUE23)
    
    # Fourth CUE Definition
    
    FinalCUE4 = (BiomassProd)/(BiomassProd - CUE22)
    
    FinalCUETab = pd.DataFrame(np.array([[FinalCUE1],[FinalCUE2],[FinalCUE3],[FinalCUE4]]), columns = ['CUE'],index = ['Def 1','Def 2','Def 3','Def 4'])
    
    
    return FinalCUETab