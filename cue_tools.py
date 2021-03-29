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