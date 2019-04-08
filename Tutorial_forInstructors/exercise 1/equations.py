#!/usr/bin/python
########################################################
# equations.py
# Author: Veronica Llorens-Rico
# Version: 1.1
# Date: December 2015
# Description: Defines differential equations of the
#   central carbon metabolism of E. coli
# Reference: Chassagnole et al, 2002
########################################################


# DIF EQUATION SYSTEM
# =======================================================
def eqs(init, t, par):
    """
    Reads the parameters from a file and describes the differential equations of the model (Chassagnole et al, 2002)
    :param init: initial conditions
    :param t: time of the integration
    :param par: parameters passed to the model
    :return: returns the diff equations to be solved by the function odeint
    """

    # Assign parameters to the different variable names
    kALDOdhap = par[0]
    kALDOeq = par[1]
    kALDOfdp = par[2]
    kALDOgap = par[3]
    kALDOgapinh = par[4]
    KDAHPSe4p = par[5]
    KDAHPSpep = par[6]
    KENOeq = par[7]
    KENOpep = par[8]
    KENOpg2 = par[9]
    KG1PATatp = par[10]
    KG1PATfdp = par[11]
    KG1PATg1p = par[12]
    KG3PDHdhap = par[13]
    KG6PDHg6p = par[14]
    KG6PDHnadp = par[15]
    KG6PDHnadphg6pinh = par[16]
    KG6PDHnadphnadpinh = par[17]
    KGAPDHeq = par[18]
    KGAPDHgap = par[19]
    KGAPDHnad = par[20]
    KGAPDHnadh = par[21]
    KGAPDHpgp = par[22]
    KPDHpyr = par[23]
    KpepCxylasefdp = par[24]
    KpepCxylasepep = par[25]
    KPFKadpa = par[26]
    KPFKadpb = par[27]
    KPFKadpc = par[28]
    KPFKampa = par[29]
    KPFKampb = par[30]
    KPFKatps = par[31]
    KPFKf6ps = par[32]
    KPFKpep = par[33]
    KPGDHatpinh = par[34]
    KPGDHnadp = par[35]
    KPGDHnadphinh = par[36]
    KPGDHpg = par[37]
    KPGIeq = par[38]
    KPGIf6p = par[39]
    KPGIf6ppginh = par[40]
    KPGIg6p = par[41]
    KPGIg6ppginh = par[42]
    KPGKadp = par[43]
    KPGKatp = par[44]
    KPGKeq = par[45]
    KPGKpg3 = par[46]
    KPGKpgp = par[47]
    KPGluMueq = par[48]
    KPGluMupg2 = par[49]
    KPGluMupg3 = par[50]
    KPGMeq = par[51]
    KPGMg1p = par[52]
    KPGMg6p = par[53]
    KPKadp = par[54]
    KPKamp = par[55]
    KPKatp = par[56]
    KPKfdp = par[57]
    KPKpep = par[58]
    KPTSa1 = par[59]
    KPTSa2 = par[60]
    KPTSa3 = par[61]
    KPTSg6p = par[62]
    KR5PIeq = par[63]
    KRPPKrib5p = par[64]
    KRu5Peq = par[65]
    KSerSynthpg3 = par[66]
    KSynth1pep = par[67]
    KSynth2pyr = par[68]
    KTAeq = par[69]
    kTISdhap = par[70]
    kTISeq = par[71]
    kTISgap = par[72]
    KTKaeq = par[73]
    KTKbeq = par[74]
    LPFK = par[75]
    LPK = par[76]
    nDAHPSe4p = par[77]
    nDAHPSpep = par[78]
    nG1PATfdp = par[79]
    nPDH = par[80]
    npepCxylasefdp = par[81]
    nPFK = par[82]
    nPK = par[83]
    nPTSg6p = par[84]
    rmaxALDO = par[85]
    rmaxDAHPS = par[86]
    rmaxENO = par[87]
    rmaxG1PAT = par[88]
    rmaxG3PDH = par[89]
    rmaxG6PDH = par[90]
    rmaxGAPDH = par[91]
    rmaxMetSynth = par[92]
    rmaxMurSynth = par[93]
    rmaxPDH = par[94]
    rmaxpepCxylase = par[95]
    rmaxPFK = par[96]
    rmaxPGDH = par[97]
    rmaxPGI = par[98]
    rmaxPGK = par[99]
    rmaxPGluMu = par[100]
    rmaxPGM = par[101]
    rmaxPK = par[102]
    rmaxPTS = par[103]
    rmaxR5PI = par[104]
    rmaxRPPK = par[105]
    rmaxRu5P = par[106]
    rmaxSerSynth = par[107]
    rmaxSynth1 = par[108]
    rmaxSynth2 = par[109]
    rmaxTA = par[110]
    rmaxTIS = par[111]
    rmaxTKa = par[112]
    rmaxTKb = par[113]
    rmaxTrpSynth = par[114]
    VALDOblf = par[115]

    # fixed parameters
    cfeed = 111.1
    Dil = 2.78e-05
    mu = 2.78e-05
    cytosol = 1
    extracellular = 1

    # Assign initial conditions
    cdhap = init[0]
    ce4p = init[1]
    cpg2 = init[2]
    cpg3 = init[3]
    cpgp = init[4]
    crib5p = init[5]
    cribu5p = init[6]
    csed7p = init[7]
    cxyl5p = init[8]
    cf6p = init[9]
    cfdp = init[10]
    cg1p = init[11]
    cg6p = init[12]
    cgap = init[13]
    cglcex = init[14]
    cpep = init[15]
    cpg = init[16]
    cpyr = init[17]

    # metabolites such as ATP, ADP, AMP... are defined by an analytical expression
    # The perturbation occurs at time=0, so if t>0 the metabolites concentration will change according to these
    # functions
    if t > 0.0:
        cadp = 0.582 + 1.73 * 2.731 ** (-0.15 * t) * (0.12 * t + 0.000214 * t ** 3)
        camp = 0.123 + 7.25 * (t / (7.25 + 1.47 * t + 0.17 * t ** 2)) + 1.073 / (1.29 + 8.05 * t)
        catp = 4.27 - 4.163 * (t / (0.657 + 1.43 * t + 0.0364 * t ** 2))
        cnad = 1.314 + 1.314 * 2.73 ** (-0.0435 * t - 0.342) - (t + 7.871) * (
               2.73 ** (-0.0218 * t - 0.171) / (8.481 + t))
        cnadh = 0.0934 + 0.00111 * 2.371 ** (-0.123 * t) * (0.844 * t + 0.104 * t ** 3)
        cnadp = 0.159 - 0.00554 * (t / (2.8 - 0.271 * t + 0.01 * t ** 2)) + 0.182 / (4.82 + 0.526 * t)
        cnadph = 0.062 + 0.332 * 2.718 ** (-0.464 * t) * (
            0.0166 * t ** 1.58 + 0.000166 * t ** 4.73 + 0.1312 * 10 ** (-9) * t ** 7.89 + 0.1362 * 10 ** (
                -12) * t ** 11 + 0.1233 * 10 ** (-15) * t ** 14.2)
    else:
        cadp = 0.582
        camp = 0.123
        catp = 4.27
        cnad = 1.314
        cnadh = 0.0934
        cnadp = 0.159
        cnadph = 0.062

    vALDO = cytosol * rmaxALDO * (cfdp - cgap * cdhap / kALDOeq) / (
        kALDOfdp + cfdp + kALDOgap * cdhap / (kALDOeq * VALDOblf) + kALDOdhap * cgap / (
            kALDOeq * VALDOblf) + cfdp * cgap / kALDOgapinh + cgap * cdhap / (VALDOblf * kALDOeq))
    vDAHPS = cytosol * rmaxDAHPS * ce4p ** nDAHPSe4p * cpep ** nDAHPSpep / (
        (KDAHPSe4p + ce4p ** nDAHPSe4p) * (KDAHPSpep + cpep ** nDAHPSpep))
    vDHAP = cytosol * mu * cdhap
    vE4P = cytosol * mu * ce4p
    vENO = cytosol * rmaxENO * (cpg2 - cpep / KENOeq) / (KENOpg2 * (1 + cpep / KENOpep) + cpg2)
    vEXTER = extracellular * Dil * (cfeed - cglcex)
    vG1PAT = cytosol * rmaxG1PAT * cg1p * catp * (1 + (cfdp / KG1PATfdp) ** nG1PATfdp) / (
        (KG1PATatp + catp) * (KG1PATg1p + cg1p))
    vG3PDH = cytosol * rmaxG3PDH * cdhap / (KG3PDHdhap + cdhap)
    vG6P = cytosol * mu * cg6p
    vG6PDH = cytosol * rmaxG6PDH * cg6p * cnadp / (
        (cg6p + KG6PDHg6p) * (1 + cnadph / KG6PDHnadphg6pinh) * (
            KG6PDHnadp * (1 + cnadph / KG6PDHnadphnadpinh) + cnadp))
    vGAP = cytosol * mu * cgap
    vGAPDH = cytosol * rmaxGAPDH * (cgap * cnad - cpgp * cnadh / KGAPDHeq) / (
        (KGAPDHgap * (1 + cpgp / KGAPDHpgp) + cgap) * (KGAPDHnad * (1 + cnadh / KGAPDHnadh) + cnad))
    vGLP = cytosol * mu * cg1p
    vMURSyNTH = cytosol * rmaxMurSynth
    vMethSynth = cytosol * rmaxMetSynth
    vPDH = cytosol * rmaxPDH * cpyr ** nPDH / (KPDHpyr + cpyr ** nPDH)
    vPEP = cytosol * mu * cpep
    vPFK = cytosol * rmaxPFK * catp * cf6p / ((catp + KPFKatps * (1 + cadp / KPFKadpc)) * (
        cf6p + KPFKf6ps * (1 + cpep / KPFKpep + cadp / KPFKadpb + camp / KPFKampb) / (
            1 + cadp / KPFKadpa + camp / KPFKampa)) * (
                                                  1 + LPFK / (1 + cf6p * (1 + cadp / KPFKadpa + camp / KPFKampa) / (
                                                      KPFKf6ps * (
                                                          1 + cpep / KPFKpep + cadp / KPFKadpb + camp / KPFKampb))) ** nPFK))
    vPG = cytosol * mu * cpg
    vPG3 = cytosol * mu * cpg3
    vPGDH = cytosol * rmaxPGDH * cpg * cnadp / (
        (cpg + KPGDHpg) * (cnadp + KPGDHnadp * (1 + cnadph / KPGDHnadphinh) * (1 + catp / KPGDHatpinh)))
    vPGI = cytosol * rmaxPGI * (cg6p - cf6p / KPGIeq) / (
        KPGIg6p * (1 + cf6p / (KPGIf6p * (1 + cpg / KPGIf6ppginh)) + cpg / KPGIg6ppginh) + cg6p)
    vPGK = cytosol * rmaxPGK * (cadp * cpgp - catp * cpg3 / KPGKeq) / (
        (KPGKadp * (1 + catp / KPGKatp) + cadp) * (KPGKpgp * (1 + cpg3 / KPGKpg3) + cpgp))
    vPGM = cytosol * rmaxPGM * (cg6p - cg1p / KPGMeq) / (KPGMg6p * (1 + cg1p / KPGMg1p) + cg6p)
    vPGP = cytosol * mu * cpgp
    vPK = cytosol * rmaxPK * cpep * (cpep / KPKpep + 1) ** (nPK - 1) * cadp / (
        KPKpep * (
            LPK * ((1 + catp / KPKatp) / (cfdp / KPKfdp + camp / KPKamp + 1)) ** nPK + (cpep / KPKpep + 1) ** nPK) * (
            cadp + KPKadp))
    vRPPK = cytosol * rmaxRPPK * crib5p / (KRPPKrib5p + crib5p)
    vPTS = extracellular * rmaxPTS * cglcex * (cpep / cpyr) / (
        (KPTSa1 + KPTSa2 * (cpep / cpyr) + KPTSa3 * cglcex + cglcex * (cpep / cpyr)) * (1 + cg6p ** nPTSg6p / KPTSg6p))
    vR5PI = cytosol * rmaxR5PI * (cribu5p - crib5p / KR5PIeq)
    vRIB5P = cytosol * mu * crib5p
    vRibu5p = cytosol * mu * cribu5p
    vRu5P = cytosol * rmaxRu5P * (cribu5p - cxyl5p / KRu5Peq)
    vSED7P = cytosol * mu * csed7p
    vSynth1 = cytosol * rmaxSynth1 * cpep / (KSynth1pep + cpep)
    vSynth2 = cytosol * rmaxSynth2 * cpyr / (KSynth2pyr + cpyr)
    vTA = cytosol * rmaxTA * (cgap * csed7p - ce4p * cf6p / KTAeq)
    vTIS = cytosol * rmaxTIS * (cdhap - cgap / kTISeq) / (kTISdhap * (1 + cgap / kTISgap) + cdhap)
    vTKA = cytosol * rmaxTKa * (crib5p * cxyl5p - csed7p * cgap / KTKaeq)
    vTKB = cytosol * rmaxTKb * (cxyl5p * ce4p - cf6p * cgap / KTKbeq)
    vTRPSYNTH = cytosol * rmaxTrpSynth
    vXYL5P = cytosol * mu * cxyl5p
    vf6P = cytosol * mu * cf6p
    vfdP = cytosol * mu * cfdp
    vpepCxylase = cytosol * rmaxpepCxylase * cpep * (1 + (cfdp / KpepCxylasefdp) ** npepCxylasefdp) / (
        KpepCxylasepep + cpep)
    vpg2 = cytosol * mu * cpg2
    vpyr = cytosol * mu * cpyr
    vrpGluMu = cytosol * rmaxPGluMu * (cpg3 - cpg2 / KPGluMueq) / (KPGluMupg3 * (1 + cpg2 / KPGluMupg2) + cpg3)
    vsersynth = cytosol * rmaxSerSynth * cpg3 / (KSerSynthpg3 + cpg3)

    dcdhap = (vALDO - vDHAP - vG3PDH - vTIS) / cytosol
    dce4p = (- vDAHPS - vE4P + vTA - vTKB) / cytosol
    dcf6p = (- 2.0 * vMURSyNTH - vPFK + vPGI + vTA + vTKB - vf6P) / cytosol
    dcfdp = (- vALDO + vPFK - vfdP) / cytosol
    dcg1p = (- vG1PAT - vGLP + vPGM) / cytosol
    dcg6p = (- vG6P - vG6PDH - vPGI - vPGM + 64.82759 * vPTS) / cytosol
    dcgap = (vALDO - vGAP - vGAPDH - vTA + vTIS + vTKA + vTKB + vTRPSYNTH) / cytosol
    dcpep = (- vDAHPS + vENO - vPEP - vPK - 64.82759 * vPTS - vSynth1 - vpepCxylase) / cytosol
    dcpg = (vG6PDH - vPG - vPGDH) / cytosol
    dcpg2 = (- vENO - vpg2 + vrpGluMu) / cytosol
    dcpg3 = (- vPG3 + vPGK - vrpGluMu - vsersynth) / cytosol
    dcpgp = (vGAPDH - vPGK - vPGP) / cytosol
    dcpyr = (vMethSynth - vPDH + vPK + 64.82759 * vPTS - vSynth2 + vTRPSYNTH - vpyr) / cytosol
    dcrib5p = (- vRPPK + vR5PI - vRIB5P - vTKA) / cytosol
    dcribu5p = (vPGDH - vR5PI - vRibu5p - vRu5P) / cytosol
    dcsed7p = (- vSED7P - vTA + vTKA) / cytosol
    dcxyl5p = (vRu5P - vTKA - vTKB - vXYL5P) / cytosol
    dcglcex = (vEXTER - vPTS) / extracellular

    return [dcdhap, dce4p, dcpg2, dcpg3, dcpgp, dcrib5p, dcribu5p, dcsed7p, dcxyl5p, dcf6p, dcfdp, dcg1p, dcg6p,
            dcgap, dcglcex, dcpep, dcpg, dcpyr]
    # =======================================================
