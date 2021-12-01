import gurobipy as gp
import time
import numpy as np

"""
Created on 10/2021

@author: Tim Adamek
"""

#berechnet das Netzwerk und verteilt die Verluste je nach Allokationsart auf

def Netzberechnung(Netzdaten, nodes, options, week, net_data, params, z, Ubergabedaten, counter):
    now = time.time()
    time_hor = options["time_hor"]
    time_steps = list(range(0, options["time_hor"]))

    mod = gp.Model()



    power = {}
    power["from_grid"] = {}  # Residual network electricity demand
    power["to_grid"] = {}  # Residual feed in
    for t in time_steps:
        power["from_grid"][t] = {}
        power["to_grid"][t] = {}

    residual = {}
    residual["power"] = {}  # Residual network electricity demand
    residual["feed"] = {}  # Residual feed in
    for t in time_steps:
            residual["power"][t] = {}
            residual["feed"][t] = {}


    # Übergabe der Daten aus der verteilten Optimierung. hier: Interaktionen mit dem übergeordneten Netz
    for t in time_steps:
        power["from_grid"][t] = sum(Netzdaten["P_from_Grid"][n][t] for n in nodes)
        power["to_grid"][t] = sum(Netzdaten["P_to_Grid"][n][t] for n in nodes)


    power_net = {}
    power_net["Inj"] = {}
    power_net["Load"] = {}
    for node in net_data["gridnodes"]:
         # Electricity demand
            power_net["Inj"][node] = {}
            power_net["Load"][node] = {}
            for t in time_steps:
                power_net["Inj"][node][t] = mod.addVar(vtype="C", name="powerInj_" + str(z) + "_t" + str(t))
                power_net["Load"][node][t] = mod.addVar(vtype="C", name="powerLoad_" + str(z) + "_t" + str(t))

    # Netzwerkleistung für einen Abschnitt [n, m]
    powerLine = mod.addVars(net_data["nodeLines"], time_steps, vtype="C", lb=-10000, name="powerLine_")
    powerLineSS = mod.addVars(net_data["nodeLines"], time_steps, vtype="C", lb=-10000, name="powerLine_")
    # set trafo bounds due to technichal limits
    powerTrafoLoad = mod.addVars(time_steps, vtype="C", lb=0, ub=net_data["trafo_max"], name="powerTrafo_" + str(t))
    powerTrafoInj = mod.addVars(time_steps, vtype="C", lb=0, ub=net_data["trafo_max"], name="injTrafo_" + str(t))

    obj = mod.addVar(vtype="C", lb=-gp.GRB.INFINITY, name="obj")
    posten = mod.addVar(vtype="C", lb=-gp.GRB.INFINITY, name="posten")

    ##Neu für Netzverluste
    yTrafo = mod.addVars(time_steps, vtype="B", name="yTrafo_" + str(t))
    power = {}
    for device in ["from_grid", "to_grid"]:
        power[device] = {}
        for t in time_steps:
            power[device][t] = mod.addVar(vtype="C", lb=0, name="power_" + device + "_t" + str(t))

    residual = {}
    residual["power"] = {}  # Residual network electricity demand
    residual["feed"] = {}  # Residual feed in
    for t in time_steps:
        residual["power"][t] = mod.addVar(vtype="C", name="residual_power_t" + str(t))
        residual["feed"][t] = mod.addVar(vtype="C", name="residual_feed_t" + str(t))

    P_loss = mod.addVars(net_data["nodeLines"], time_steps, vtype="C", lb=-10000, name="powerLoss")
    P_loss_Ges = mod.addVar(vtype="C", lb=-gp.GRB.INFINITY, name="obj")

    I_quad = mod.addVars(net_data["nodeLines"], time_steps, vtype="C", lb=-10000, name="Current_quad")
    Res = mod.addVars(net_data["nodeLines"], vtype="C", lb=-10000, name="Resistance")
    Resqua = mod.addVars(net_data["nodeLines"], vtype="C", lb=-10000, name="Resistancequad")
    V_HIlF = mod.addVars(net_data["nodeLines"], time_steps, vtype="C", lb=-10000, name="VHilfe")
    V_quad = {}
    for i in range(1, 228):
        V_quad[i] = {}
        for t in time_steps:
            V_quad[i][t] = mod.addVar(vtype="C", name="VoltageQuad" + str(i) + "_t" + str(t))



    mod.update()
    mod.setObjective(posten)
    mod.ModelSense = gp.GRB.MINIMIZE

    # %% ENERGY BALANCES: NETWORK AND ENERGY HUB

    # for t in time_steps:
    # For all modes and scenarios:
    # Electricity balance (Power balance, MW)




    # %% NETWORK CONSTRAINTS
    for t in time_steps:
        mod.addConstr(residual["feed"][t] == sum(Netzdaten["P_inj_Ges"][n][t] for n in nodes))
        mod.addConstr(residual["power"][t] == sum(Netzdaten["P_load_Ges"][n][t] for n in nodes))

        # Energiebilanz für gesamtes Netz
        mod.addConstr(residual["feed"][t] + power["from_grid"][t] - sum(P_loss[n, m, t] for [n, m] in net_data["nodeLines"]) ==residual["power"][t] + power["to_grid"][t], name="NBEBGesamt")

        mod.addConstr(power["to_grid"][t] <= residual["feed"][t])

        mod.addConstr(power["from_grid"][t] <= yTrafo[t] * 1000, name="Binary1_"  + "_" + str(t))
        mod.addConstr(power["to_grid"][t] <= (1 - yTrafo[t]) * 1000, name="Binary2_" + "_" + str(t))

    for [n, m] in net_data["nodeLines"]:
        if [n, m] == [1, 0]:
            pass
        else:

            mod.addConstr(Res[n, m] == net_data["res"][n, m] * net_data["length"][n, m], name="1") # Berechnung Widerstand = spez. Widerstand * Länge
            mod.addConstr(Resqua[n, m] == net_data["resQua"][n, m], name="2") # quadratischer Widerstand
            #Res[n, m]**2
            # LinDistFlow Gleichungen ohne Nebenbedingung zur Spannungshaltung -> wird am Ende überprüft, da hier bereits feste Werte für Handelsmengen
            for t in time_steps:
                mod.addConstr(P_loss[n, m, t] == Res[n, m] * I_quad[n, m, t] / 1000, name="3")

                mod.addConstr(V_quad[m][t] == V_quad[n][t] - 2 * Res[n, m] * powerLine[n, m, t] / 1000 + V_HIlF[
                    n, m, t] / 1000000, name="4")
                mod.addConstr(V_HIlF[n, m, t] == ((Resqua[n, m]) * I_quad[n, m, t]), name="5")
                mod.addGenConstrAbs(powerLineSS[n, m, t], powerLine[n, m, t], "absconstr")

                # konvexe Relaxation
                mod.addConstr((0.4)**2 * I_quad[n, m, t] >= (powerLine[n, m, t]) ** 2, name="6")
                # V_quad[n][t]


    for t in time_steps:
        mod.addConstr(V_quad[1][t] == (0.400) ** 2, name="9")

    #for node in net_data["gridnodes"]:
        #if node == 0:
            #pass
        #else:
            #for t in time_steps:
                #mod.addConstr(V_quad[node][t] <= (0.440) ** 2, name="7")
                #mod.addConstr(V_quad[node][t] >= (0.360) ** 2, name="8")


    #Gesamtverluste
    mod.addConstr(P_loss_Ges == sum(P_loss[n, m, t] for [n, m] in net_data["nodeLines"] for t in time_steps))






    for t in time_steps:
            powerTrafoLoad[t] == power["from_grid"][t]
            powerTrafoInj[t] == power["to_grid"][t]

    for t in time_steps:

            for index, row in net_data["grid_allo"].iterrows():
                mod.addConstr(power_net["Inj"][row["gridnodes"]][t] == Netzdaten["P_inj_Ges"][row["nodes"]][t])
                mod.addConstr(power_net["Load"][row["gridnodes"]][t] == Netzdaten["P_load_Ges"][row["nodes"]][t])

    for node in net_data["gridnodes"]:

             for t in time_steps:

                if node in net_data["net_nodes"]["load"]:

                    dummy = 1

                else:

                    mod.addConstr(power_net["Inj"][node][t] == 0)
                    mod.addConstr(power_net["Load"][node][t] == 0)

    for node in net_data["gridnodes"]:

             for t in time_steps:

                if node in net_data["net_nodes"]["trafo"]:
                    # Knotengleichung mit Trafobeteiligung
                    mod.addConstr(powerLine.sum(node, '*', t) - powerLine.sum('*', node, t) ==
                                powerTrafoLoad[t] - powerTrafoInj[t], name="node balance_" + str(z))

                    mod.addConstr(power_net["Inj"][node][t] == 0)
                    mod.addConstr(power_net["Load"][node][t] == 0)

                else:
                    # Knotengleichung
                    mod.addConstr(powerLine.sum(node, '*', t) - powerLine.sum('*', node, t) ==
                                power_net["Inj"][node][t] - power_net["Load"][node][t] - P_loss.sum('*',node,t),
                                name="node balance_" + str(node))

    for [n, m] in net_data["nodeLines"]:

                for t in time_steps:
                    mod.addConstr(powerLine[n, m, t] <= net_data["powerLine_max"][n, m],
                            name="line power max_" + str(z) + str(m) + str(t))
                    mod.addConstr(powerLine[n, m, t] >= (-1) * net_data["powerLine_max"][n, m],
                            name="line power min_" + str(z) + str(m) + str(t))

    #  Zielfunktion Gesamtverluste
    mod.addConstr(posten == P_loss_Ges)


    mod.addConstr(obj == posten)

    #%% START OPTIMIZATION
    # set objective function
    mod.setObjective(obj, gp.GRB.MINIMIZE)

    mod.Params.NonConvex = 2
    mod.Params.NumericFocus = 3
    mod.Params.MIPGap = 0.08

    mod.optimize()

    if mod.status == gp.GRB.Status.INFEASIBLE:
        print('is infisible')

    later = time.time()
    # Berechnung Rechenzeit Netzwerk
    DeltaT = later - now


    print("********************************************")
    print("GridModel run time was " + str(DeltaT) + " seconds")
    print("********************************************")



    print(P_loss_Ges.X, 'Verlust')
        #for t in time_steps:
            #if powerLine[n, m, t] >= net_data["powerLine_max"][n, m]:
                    #and powerLine[n, m, t] <= (-1) * net_data["powerLine_max"][n, m]:
                #print('FEHLER')
            #model.addConstr(powerLine[n, m, t] <= net_data["powerLine_max"][n, m],
                         #   name="line power max_" + str(z) + str(m) + str(t))
            #model.addConstr(powerLine[n, m, t] >= (-1) * net_data["powerLine_max"][n, m],
                           # name="line power min_" + str(z) + str(m) + str(t))



    # Allocation due to power
    if options["Allocation"]:
        ule = sum(Netzdaten["P_P2P_inj"][n][j][t]  for n in nodes for j in nodes for t in time_steps) + sum(Netzdaten["P_to_Grid"][n][t] for n in nodes for t in time_steps ) + sum(Netzdaten["P_from_Grid"][n][t] for n in nodes for t in time_steps )
        Verl ={}
        for n in nodes:
            Verl[n] = {}
            for j in nodes:
                Verl[n][j] = {}
                for t in time_steps:
                    Verl[n][j][t] = P_loss_Ges.X * Netzdaten["P_P2P_inj"][n][j][t]/(ule)
    else:
    # Allocation due to power^2*Resistance

        kk = sum(Netzdaten["P_P2P_inj"][n][j][t]**2 * net_data["om"][n][j] for n in nodes for j in nodes for t in time_steps) + sum((Netzdaten["P_to_Grid"][n][t])**2 * net_data["omGrid"][n] for n in nodes for t in time_steps) + sum((Netzdaten["P_from_Grid"][n][t])**2 * net_data["omGrid"][n] for n in nodes for t in time_steps)
        Verl ={}
        for n in nodes:
            Verl[n] = {}
            for j in nodes:
                Verl[n][j] = {}
                for t in time_steps:
                    Verl[n][j][t] = P_loss_Ges.X * Netzdaten["P_P2P_inj"][n][j][t]**2 * net_data["om"][n][j]/(kk)

    Wahrheit = True

    return Wahrheit, Verl, DeltaT