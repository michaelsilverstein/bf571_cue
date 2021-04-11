import tempfile
import cobra
import requests
import os
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
    
    ## knock-outs loop ##
    # initialize empty matricies
    rxns = len(model.reactions)
    
    v1_grow = np.zeros((rxns,rxns))
    v2_grow = np.zeros((rxns,rxns))
    v1v2_grow = np.zeros((rxns,rxns))
    for i in range(rxns):
        for j in range(rxns):
            if j > i:
                # Buffer variables
                upper_i = model.reactions[i].upper_bound
                lower_i = model.reactions[i].lower_bound
                upper_j = model.reactions[j].upper_bound
                lower_j = model.reactions[j].lower_bound

                ## first reaction knock out ##
                # set upper and lower bounds to zero
                model.reactions[i].upper_bound = 0
                model.reactions[i].lower_bound = 0
                # solve model and record growth rate
                solution = model.optimize()
                v1_grow[i,j] = solution.objective_value
                # return bounds to their previous state
                model.reactions[i].upper_bound = upper_i
                model.reactions[i].lower_bound = lower_i

                ## Second rection knock out ##
                model.reactions[j].upper_bound = 0
                model.reactions[j].lower_bound = 0
                # Solve model and record growth rate
                solution = model.optimize()
                v2_grow[i,j] = solution.objective_value

                ## combo knock out ##
                # Set bounds on other rxn to zero, now both are zero
                model.reactions[i].upper_bound = 0
                model.reactions[i].lower_bound = 0
                # Solve
                solution = model.optimize()
                v1v2_grow[i,j] = solution.objective_value
                # return them to what they were
                model.reactions[i].upper_bound = upper_i
                model.reactions[i].lower_bound = lower_i
                model.reactions[j].upper_bound = upper_j
                model.reactions[j].lower_bound = lower_j

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

