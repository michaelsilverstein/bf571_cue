with model:
    medium = model.medium
    medium["EX_co2_e"] = 1000
    medium["EX_glc__D_e"] = 10
    #medium["EX_glu__L_e"] = 0000
    #medium["EX_gln__L_e"] = 1000
    #medium["EX_fru_e"] = 1000
    #medium["EX_mal__L_e"] = 0000
    #medium["EX_ac_e"] = 10
    medium["EX_h2o_e"] = 1000
    medium["EX_h_e"] = 1000
    medium["EX_nh4_e"] = 4
    medium["EX_o2_e"] = 1000
    medium["EX_pi_e"] = 1000
    #medium[""]
    model.medium = medium
    #print(model.optimize())
    #print('%f', model.reactions.BIOMASS_Ecoli_core_w_GAM)
    #cue = computeCUE(model, biomass_rxn, co2_secretion)
    cue = CUEDef(model)
    print(cue)
    print("%s ; %f" % (model.reactions.PGI.id, model.reactions.PGI.flux))
    print(model.summary())
