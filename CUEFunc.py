#!/usr/bin/env python
# coding: utf-8

# # CUE Function 
# ### This function will take a cobrapy model and give you the CUE based on four definitions
# 

# In[7]:


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


# In[ ]:




