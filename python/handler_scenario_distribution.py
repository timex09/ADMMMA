#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 10 2021

@author: Tim Adamek
"""
# import python packages
import numpy as np
from numpy.linalg import norm

# import own functions
import python.opti_stochastic_distributed as stoch_opti_dist
import python.GridBerucksichtigung as BerechnungNetz


# ADMM verteilte Optimierung für Energy-Communites

def admm(nodes, options, week, net_data, params, result_file):




    # Time Parameter
    time_hor = options["time_hor"]
    rounds = {}
    Speicher={}
    costs = {}
    res = {}
    zeit={}
    optiparams = {}


    gesamtzeit=np.zeros(4) # Erstellung Gesamtzeit als Nullvektor mit 4 Einträgen für 4 Typtage

    time_steps = list(range(0, options["time_hor"]))
    dt = 1
    last_time_step = options["time_hor"]
    idx = week


    #%% Optimization Iteration
    status = False    # Status der Konvergenz, am Anfang = 0
    counter = 0     # Iteration counter am Anfang 0 = 1. Iterationsschritt

    AA=np.zeros((600, len(time_steps)))
    HilfsvariableA = AA
    HilfsvariableB = AA
    HilfsvariableC = AA
    HilfsvariableD = AA



    primal=np.zeros((len(nodes), len(time_steps)))
    dual=np.zeros((len(nodes), len(time_steps)))


    # penalty parameter rho: Startwertübergabe
    rho = np.zeros(600)
    for i in range(0, 600):
        rho[i]=options["rho"]

    # Preis Lambda
    lam={}
    for n in nodes:
        lam[n]={}
        for t in time_steps:
            lam[n][t]={}

    # hilfsvariablen
    PP=np.zeros(( len(nodes) , len(time_steps) ))
    FF=np.zeros(( len(nodes) , len(time_steps) ))
    WW=np.zeros(( len(nodes) , len(time_steps) ))
    EE=np.zeros(( len(nodes) , len(time_steps) ))
    OO=np.zeros(( len(nodes) , len(time_steps) ))
    II=np.zeros(( len(nodes) , len(time_steps) ))
    LKA=np.zeros(( len(nodes) ,len(nodes),  len(time_steps) ))
    PVLeistung = np.zeros((len(nodes), len(time_steps)))

    JJ=np.zeros((4 , 700 ))


    MaxZeit=np.zeros((len(nodes)))
    maximumZeit=np.zeros((4 , 700 ))
    NetzberechZeit = np.zeros((4 , 700 ))

    ZfWert=np.zeros((4 , 700 ))



    Speicher["p_load_P2P"] = np.zeros((600, len(nodes), len(nodes) , len(time_steps) ))
    Speicher["p_inj_P2P"] = np.zeros((600, len(nodes), len(nodes) , len(time_steps) ))

    P_P2P_inj=np.zeros((len(nodes), len(nodes) , len(time_steps) ))
    P_P2P_load=np.zeros((len(nodes), len(nodes) , len(time_steps) ))


    HilfsvariableE = np.zeros((600, len(nodes), len(time_steps)))
    HilfsvariableF = np.zeros((600, len(nodes), len(time_steps)))


    # Erstellung der Übergabedaten
    Ubergabedaten = {"P_P2P_inj":np.zeros((len(nodes), len(nodes) , len(time_steps) )),
                         "P_P2P_load":np.zeros((len(nodes), len(nodes) , len(time_steps) )),
                         "P_Lade_Netz":HilfsvariableA,
                         "P_Einsp_Netz":HilfsvariableB,
                         "P_Nutz":HilfsvariableC,
                         "P_Verkauf":HilfsvariableD,
                        "P_Sammel_inj":HilfsvariableE,
                        "P_Sammel_load":HilfsvariableF,
                        "P_load_netz": np.zeros(( len(nodes) , len(time_steps) )),
                        "P_inj_netz": np.zeros(( len(nodes) , len(time_steps) )),
                         "Verl": np.zeros(( len(nodes) ,len(nodes),  len(time_steps) )),

                     }

    #Erstellung der Netzdaten, später für die Netzwerkberechnung
    Netzdaten = {"P_to_Grid": np.zeros(( len(nodes) , len(time_steps) )),
                             "P_from_Grid": np.zeros(( len(nodes) , len(time_steps) )),
                             "P_inj_Ges": np.zeros(( len(nodes) , len(time_steps) )),
                             "P_load_Ges": np.zeros(( len(nodes) , len(time_steps) )),
                             "P_P2P_inj": np.zeros(( len(nodes) ,len(nodes),  len(time_steps) )),
                             }



    while status == False and counter < options["admm_iteration_limit"]:

        # 1. Iterationsschritt
        if counter == 0:
            # Startwert für die Preise Lambda
            for n in nodes:
                for t in time_steps:
                    lam[n][t]=options["admm_startval_dual"]


            # Übergabe der Preise für die verteilte Optimierung
            admm_params = { "lambda_load": lam}

            nn = 0
            # Optimierung für jeden Haushalt z
            for z in nodes:
                # startet verteilte optimierung
                costs[z], res[z], optiparams[z], zeit[z] = stoch_opti_dist.opti(nodes, options, week, net_data, params, z, Ubergabedaten, counter,admm_params, rho, Netzdaten)
                # Berechnung der Gesamtrechenzeit
                gesamtzeit[week]=gesamtzeit[week] + zeit[z]
                # Berechnung der maximalen Rechenzeit eines Haushaltes -> Optimierung gleichzeitig
                MaxZeit[z]=zeit[z]
                # Ergebnisse speichern
                rounds[counter] = res.copy()

                # Berechnen der Gesamtkosten durch Summation aller Einzelkosten
                nn = nn + costs[z]
                if z == 112:
                    print('Gesamtkosten=', nn)
                #Übergabe upgedateter Werte
                for t in time_steps:
                    HilfsvariableA[counter][t] = HilfsvariableA[counter][t] + (res[z]["p_load_grid"][t])
                    HilfsvariableB[counter][t] = HilfsvariableB[counter][t] + (res[z]["p_inj_grid"][t])
                    HilfsvariableC[counter][t] = HilfsvariableC[counter][t] + (res[z]["p_use"][t])
                    HilfsvariableD[counter][t] = HilfsvariableD[counter][t] + (res[z]["p_sell"][t])
                    HilfsvariableE[counter][z][t] = res[z]["p_sell"][t]
                    HilfsvariableF[counter][z][t] = res[z]["p_use"][t]
                    PP[z][t]=res[z]["p_load_grid"][t]
                    FF[z][t]=res[z]["p_inj_grid"][t]
                    WW[z][t] = res[z]["p_sell"][t]
                    EE[z][t] = res[z]["p_use"][t]

                    OO[z][t]=res[z]["p_inj_P2P_Ges"][t]
                    II[z][t]=res[z]["p_load_P2P_Ges"][t]

                    PVLeistung[z][t] = res[z]["PV_P"][t]

                    for j in nodes:
                        P_P2P_inj[z][j][t] = res[z]["p_inj_P2P"][j][t]
                        P_P2P_load[z][j][t] = res[z]["p_load_P2P"][j][t]
                        Speicher["p_load_P2P"][counter][z][j][t] = P_P2P_load[z][j][t]
                        Speicher["p_inj_P2P"][counter][z][j][t] = P_P2P_inj[z][j][t]




                # wenn Optimierung nacheinander hier Übergabedaten erneuern
                if options["sequence"] == False:
                    Ubergabedaten = {"P_P2P_inj":P_P2P_inj,
                         "P_P2P_load":P_P2P_load,
                         "P_Lade_Netz":HilfsvariableA,
                         "P_Einsp_Netz":HilfsvariableB,
                         "P_Nutz":HilfsvariableC,
                         "P_Verkauf":HilfsvariableD,
                        "P_Sammel_inj": HilfsvariableE,
                        "P_Sammel_load": HilfsvariableF,
                        "P_load_netz": PP,
                        "P_inj_netz": FF,
                                     }

                # Übergabedaten für Netzwerkberechnung
                Netzdaten = {"P_to_Grid": FF,
                                 "P_from_Grid": PP,
                                 "P_inj_Ges": WW,
                                 "P_load_Ges": EE,
                             "P_P2P_inj": P_P2P_inj,
                                 }

            # berechnet primales Residuum, kein duales, da 1. Iterationsschritt
            for n in nodes:
                for t in time_steps:
                        primal[n][t] = (sum((P_P2P_load[j][n][t] - P_P2P_inj[n][j][t]) for j in nodes))


            # Status falsch, da 1. Iterationsschritt
            status=False
        else:

            print("********************************")
            print("This is iteration No." + str(counter))
            print("********************************")

            nn = 0
            for z in nodes:

                costs[z], res[z], optiparams[z], zeit[z]= stoch_opti_dist.opti(nodes, options, week, net_data, params, z, Ubergabedaten, counter,admm_params, rho, Netzdaten)
                gesamtzeit[week] = gesamtzeit[week] + zeit[z]
                MaxZeit[z] = zeit[z]

                # Safe results of each admm round
                rounds[counter] = res.copy()

                nn=nn+costs[z]
                if z==112:
                    print('Gesamtkosten=', nn)
                for t in time_steps:
                    HilfsvariableA[counter][t] = HilfsvariableA[counter][t] + (res[z]["p_load_grid"][t] )
                    HilfsvariableB[counter][t] = HilfsvariableB[counter][t] + (res[z]["p_inj_grid"][t])
                    HilfsvariableC[counter][t] = HilfsvariableC[counter][t] + (res[z]["p_use"][t])
                    HilfsvariableD[counter][t] = HilfsvariableD[counter][t] + (res[z]["p_sell"][t])
                    HilfsvariableE[counter][z][t] = res[z]["p_sell"][t]
                    HilfsvariableF[counter][z][t] = res[z]["p_use"][t]
                    PP[z][t] = res[z]["p_load_grid"][t]
                    FF[z][t] = res[z]["p_inj_grid"][t]
                    WW[z][t] = res[z]["p_sell"][t]
                    EE[z][t] = res[z]["p_use"][t]

                    OO[z][t] = res[z]["p_inj_P2P_Ges"][t]
                    II[z][t] = res[z]["p_load_P2P_Ges"][t]

                    PVLeistung[z][t] = res[z]["PV_P"][t]

                    for j in nodes:
                        P_P2P_inj[z][j][t] = res[z]["p_inj_P2P"][j][t]
                        P_P2P_load[z][j][t] = res[z]["p_load_P2P"][j][t]
                        Speicher["p_load_P2P"][counter][z][j][t] = P_P2P_load[z][j][t]
                        Speicher["p_inj_P2P"][counter][z][j][t] = P_P2P_inj[z][j][t]

                # wenn Optimierung nacheinander hier Übergabedaten erneuern
                if options["sequence"] == False:
                    Ubergabedaten = {"P_P2P_inj": P_P2P_inj,
                                     "P_P2P_load": P_P2P_load,
                                     "P_Lade_Netz": HilfsvariableA,
                                     "P_Einsp_Netz": HilfsvariableB,
                                     "P_Nutz": HilfsvariableC,
                                     "P_Verkauf": HilfsvariableD,
                                     "P_Sammel_inj": HilfsvariableE,
                                     "P_Sammel_load": HilfsvariableF,
                                     "P_load_netz": PP,
                                     "P_inj_netz": FF,
                                     }
                # Übergabedaten für Netzwerkberechnung
                Netzdaten = {"P_to_Grid": FF,
                             "P_from_Grid": PP,
                             "P_inj_Ges": WW,
                             "P_load_Ges": EE,
                             "P_P2P_inj": P_P2P_inj,
                             }




            # Berechnung beider Residuen
            for n in nodes:
                for t in time_steps:
                        primal[n][t] = (sum((P_P2P_load[j][n][t] - P_P2P_inj[n][j][t]) for j in nodes))
                        dual[n][t] = abs(rho[counter] * sum((Speicher["p_load_P2P"][counter][n][j][t] - Speicher["p_load_P2P"][counter-1][n][j][t] + Speicher["p_inj_P2P"][counter][n][j][t] - Speicher["p_inj_P2P"][counter-1][n][j][t]) for j in nodes))

            # Prüfe, ob für alle Punkte Residuen eingehalten werden
            status_all = np.full((len(nodes), len(time_steps)), False, dtype=bool)
            x=0 # Anzahl richtiger Punkte
            y=0 # Anzahl falscher Punkte
            for n in nodes:
                for t in time_steps:
                    if abs(primal[n][t]) > options["admm_threshold"] or dual[n][t] > options["admm_threshold"]:

                        status_all[n][t] = False
                        x=x+1

                    else:
                        status_all[n][t] = True
                        y=y+1

            JJ[week][counter] = y

            print('falsch:', x)
            print('wahr:', y)

            # Status nur = True, wenn für alle Punkte beide Residuen eingehalten werden
            status = status_all.all()
            if counter > 140:
                status = True
            if x < 3:
                status = True


            # Berechnung Norm der Residuen
            primal_Ges = norm(primal[n][t])
            dual_Ges = norm(dual[n][t])

            # update rho
            if abs(primal_Ges) > 10*abs(dual_Ges):
                rho[counter+1]=2*rho[counter]
            elif 10*abs(primal_Ges) < abs(dual_Ges):
                rho[counter + 1] =  (rho[counter])/2
            else:
                rho[counter + 1] = rho[counter]
            # Einhalten rho < 0,5
            if rho[counter+1]<0.5:
                rho[counter + 1]=0.5






        # Update der Preise Lambda 0.008 = 1/113

        for t in time_steps:
            for n in nodes:
                lam[n][t] = lam[n][t] + 0.008*rho[counter] * (primal[n][t])

        # Netzberechnung
        if options["grid"]:
            Ja, Verl, DeltaT = BerechnungNetz.Netzberechnung(Netzdaten, nodes, options, week, net_data, params, z, Ubergabedaten, counter)

            # Übergabe der Verluste
            for n in nodes:
                for j in nodes:
                    for t in time_steps:
                        LKA[n][j][t] = Verl[n][j][t]
            # Übergabe Berechnungszeit der Netzwerkberechnung
            NetzberechZeit[week][counter] = DeltaT

        Gesamtverlust=sum(Verl[n][j][t] for n in nodes for j in nodes for t in time_steps)
        GesamtEinsp=sum(Netzdaten["P_inj_Ges"][n][t] for n in nodes for t in time_steps)
        GesamtBezug = sum(Netzdaten["P_load_Ges"][n][t] for n in nodes for t in time_steps)





        # count up

        counter = counter + 1

        # falls Optimierung der Haushalte gleichzeitg, hier Übergabedaten aktualisieren
        if options["sequence"] == True:
            Ubergabedaten = {"P_P2P_inj": P_P2P_inj,
                             "P_P2P_load": P_P2P_load,
                             "P_Lade_Netz": HilfsvariableA,
                             "P_Einsp_Netz": HilfsvariableB,
                             "P_Nutz": HilfsvariableC,
                             "P_Verkauf": HilfsvariableD,
                             "P_Sammel_inj": HilfsvariableE,
                             "P_Sammel_load": HilfsvariableF,
                             "P_load_netz": PP,
                             "P_inj_netz": FF,
                             "Verl": LKA,
                             }


        # maximale Berechnungszeit eines Haushaltes je Iterationsschrittes
        maximumZeit[week][counter] = max(MaxZeit)

        # Übergabe Gesamtkosten in jedem Iterationsschritt, um Verlauf nachvollziehen zu können
        ZfWert[week][counter] = nn

        # Übergabe Preise Lambda
        admm_params = {"lambda_load": lam}










    # Übergabe von Daten, die später als Excel ausgegeben werden
    result_file["1"] = [JJ[0][z] for z in range(700)]
    result_file["20"] = [Netzdaten["P_to_Grid"][n][t] for n in nodes for t in time_steps]
    result_file["21"] = [Netzdaten["P_from_Grid"][n][t] for n in nodes for t in time_steps]
    result_file["22"] = [II[n][t] for n in nodes for t in time_steps]
    result_file["23"] = [OO[n][t] for n in nodes for t in time_steps]
    result_file["24"] = [lam[n][t] for n in nodes for t in time_steps]
    result_file["25"] = [PVLeistung[n][t] for n in nodes for t in time_steps]
    if week==0:
        result_file["26"] = [sum(maximumZeit[week][u] for u in range(700))]
        result_file["30"] = [ZfWert[week][u] for u in range(700)]
        result_file["NetzZeit1"] = [sum(NetzberechZeit[week][u] for u in range(700))]


    elif week==1:
        result_file["27"] = [sum(maximumZeit[week][u] for u in range(700))]
        result_file["31"] = [ZfWert[week][u] for u in range(700)]
        result_file["NetzZeit2"] = [sum(NetzberechZeit[week][u] for u in range(700))]
    elif week==2:
        result_file["28"] = [sum(maximumZeit[week][u] for u in range(700))]
        result_file["32"] = [ZfWert[week][u] for u in range(700)]
        result_file["NetzZeit3"] = [sum(NetzberechZeit[week][u] for u in range(700))]
    elif week==3:
        result_file["29"] = [sum(maximumZeit[week][u] for u in range(700))]
        result_file["33"] = [ZfWert[week][u] for u in range(700)]
        result_file["NetzZeit4"] = [sum(NetzberechZeit[week][u] for u in range(700))]



    result_file["Berechnungszeit"] = gesamtzeit[week]
    result_file["Iterationsschritt"] = counter
    result_file["Gesamtkosten"] = nn
    result_file["Haushalte"] = len(nodes)




    if status==True:
        if week==0:
            result_file["5"] = [gesamtzeit[0]]
            result_file["1Verl"] = [Gesamtverlust]
        elif week==1:
            result_file["6"] = [gesamtzeit[1]]
            result_file["2Verl"] = [Gesamtverlust]
        elif week==2:
            result_file["7"] = [gesamtzeit[2]]
            result_file["3Verl"] = [Gesamtverlust]
        elif week==3:
            result_file["8"] = [gesamtzeit[3]]
            result_file["4Verl"] = [Gesamtverlust]

    if status == True:
        if week==0:
            result_file["9"] = [lam[n][t] for n in nodes for t in time_steps]
        elif week==1:
            result_file["10"] = [lam[n][t] for n in nodes for t in time_steps]
        elif week==2:
            result_file["11"] = [lam[n][t] for n in nodes for t in time_steps]
        elif week==3:
            result_file["12"] = [lam[n][t] for n in nodes for t in time_steps]

    if status == True:
        if week==0:
            result_file["13"] = [nn]
        elif week==1:
            result_file["14"] = [nn]



    return counter, status, costs, res, optiparams, rounds, result_file



