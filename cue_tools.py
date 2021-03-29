import tempfile
import cobra
import requests
import os

def loadModelURL(url):
    """
    Load SMBL cobra model from a URL
    
    Example: 
    url = 'http://bigg.ucsd.edu/static/models/e_coli_core.xml'
    model = loadModelURL(url)
    """
    with tempfile.NamedTemporaryFile('w', delete=False) as tf:
        handle = requests.get(url)
        tf.write(handle.text)
        model = cobra.io.read_sbml_model(tf.name)
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
    url = 'http://bigg.ucsd.edu/static/models/e_coli_core.xml'
    model = loadModelURL(url)
    # Calculate cue
    cue = computeCUE(model, 'BIOMASS_Ecoli_core_w_GAM', 'EX_co2_e')
    """
    # Run FBA
    res = model.optimize()
    # Compute CUE
    cue = res.get_primal_by_id(biomass_rxn)/res.get_primal_by_id(co2_rxn)
    return cue