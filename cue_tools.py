def loadModelURL(url):
    """Load SMBL cobra model from a URL"""
    with tempfile.NamedTemporaryFile('w') as tf:
        handle = requests.get(url)
        tf.write(handle.text)
        model = cobra.io.read_sbml_model(tf.name)
        return model