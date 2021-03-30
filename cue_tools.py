import tempfile
import cobra
import requests
import os

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