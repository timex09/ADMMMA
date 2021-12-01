#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 10 2021

@author: Tim Adamek
"""

import gurobipy as gp
import numpy as np
import time
import pandas as pd
import python.store_items as store

def opti(nodes, options, week, net_data, params, z,Ubergabedaten, counter,admm_params, rho, Netzdaten):
    now = time.time()

        
    #%% OPTIMIZATION MODEL

    model = gp.Model("stochastic_opti")
    model.Params.MIPGap = 0.05
    time_steps = list(range(0, options["time_hor"]))
    dt = 1
    last_time_step = options["time_hor"]
    idx = week

    # %% TECHNICAL PARAMETERS
    soc_nom = {}
    soc_nom["TES"] = {}
    soc_nom["BAT"] = {}
    soc_nom["EV"] = {}
    soc_init = {}
    soc_init["TES"] = {}
    soc_init["BAT"] = {}
    soc_init["EV"] = {}

        # Define nominal SOC_nom according to capacities: SOC: Ladezustand 0,1=10%
    soc_nom["TES"][z] = params["phy"]["c_w"] * nodes[z]["TES"]["cap"] * \
                            nodes[z]["TES"]["dT_max"] * params["phy"]["rho_w"] / 3600
    soc_nom["BAT"][z] = nodes[z]["BAT"]["cap"]
    if options["ev_status"]:
            soc_nom["EV"][z] = nodes[z]["EV"]["cap"]  # kWh   Nomincal capcity Battery from ev
    else:
            soc_nom["EV"][z] = 0
        # Initial state of charge
    soc_init["TES"][z] = soc_nom["TES"][z] * 0.5  # kWh   Initial SOC TES
    soc_init["BAT"][z] = soc_nom["BAT"][z] * 0.5  # kWh   Initial SOC Battery
    soc_init["EV"][z] = soc_nom["EV"][z] * 0.75

        # %% TECHNICAL CONSTRAINTS

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # CREATE VARIABLES

    # %% OPERATIONAL BUILDING VARIABLES
    #   NOTE: Subscript "dom" = "domestic"
    # Electrical power to/from electricity-based domestic devices
    power_dom = {}
    for device in ["HP", "EB", "PV"]:
        power_dom[device] = {}

        power_dom[device][z] = {}
        for t in time_steps:
                power_dom[device][z][t] = model.addVar(vtype="C",
                                                       name="power_" + device + "_n" + str(z) + "_t" + str(t))

    # Heat to/from devices
    heat_dom = {}
    for device in ["HP", "EB"]:
        heat_dom[device] = {}

        heat_dom[device][z] = {}
        for t in time_steps:
                heat_dom[device][z][t] = model.addVar(vtype="C", name="heat_" + device + "_n" + str(z) + "_t" + str(t))

    # Storage variables
    soc_dom = {}  # State of charge
    ch_dom = {}
    dch_dom = {}
    for device in ["TES", "BAT", "EV"]:
        soc_dom[device] = {}  # Energy (kWh)
        ch_dom[device] = {}  # Power(kW)
        dch_dom[device] = {}  # Power (kW)

        soc_dom[device][z] = {}
        ch_dom[device][z] = {}
        dch_dom[device][z] = {}
        for t in time_steps:
                soc_dom[device][z][t] = model.addVar(vtype="C", name="soc_" + device + "_n" + str(z) + "_t" + str(t))
                ch_dom[device][z][t] = model.addVar(vtype="C", name="ch_dom_" + device + "_n" + str(z) + "_t" + str(t))
                dch_dom[device][z][t] = model.addVar(vtype="C",
                                                     name="dch_dom_" + device + "_n" + str(z) + "_t" + str(t))

    # Residual building demands (in kW)  [Sum for each building of all devices]
    res_dom = {}
    res_dom["power"] = {}
    res_dom["feed"] = {}
    res_dom["power"][z] = {}
    res_dom["feed"][z] = {}

    for t in time_steps:
            res_dom["power"][z][t] = model.addVar(vtype="C", name="residual_power_n" + str(z) + "_t" + str(t))
            res_dom["feed"][z][t] = model.addVar(vtype="C", name="residual_feed_n" + str(z) + "_t" + str(t))


    if options["grid"]:
        res_dom["LoadwithLoss"] = {}
        res_dom["InjwithLoss"] = {}
        res_dom["LoadwithLoss"][z] = {}
        res_dom["InjwithLoss"][z] = {}
        for t in time_steps:
            res_dom["LoadwithLoss"][z][t] = model.addVar(vtype="C", name="residual_power_n" + str(z) + "_t" + str(t))
            res_dom["InjwithLoss"][z][t] = model.addVar(vtype="C", name="residual_feed_n" + str(z) + "_t" + str(t))

    # P2P trading quantitites: grid part
    power_grid = {}
    power_grid["inj"] = {}
    power_grid["load"] = {}
    power_grid["inj"][z] = {}
    power_grid["load"][z] = {}
    for t in time_steps:
            power_grid["inj"][z][t] = model.addVar(vtype="C", name="power_from_grid" + str(z) + "_t" + str(t))
            power_grid["load"][z][t] = model.addVar(vtype="C", name="power_to_grid" + str(z) + "_t" + str(t))

    # P2P trading quantities: P2P injection and load
    power_P2P = {}
    power_P2P["inj"] = {}
    power_P2P["load"] = {}

    power_P2P["inj"][z] = {}
    power_P2P["load"][z] = {}
    for j in nodes:
            power_P2P["inj"][z][j] = {}
            power_P2P["load"][z][j] = {}
            for t in time_steps:
                power_P2P["inj"][z][j][t] = model.addVar(vtype="C", name="power_from_P2P" + str(z) + "_t" + str(t))
                power_P2P["load"][z][j][t] = model.addVar(vtype="C", name="power_to_P2P" + str(z) + "_t" + str(t))

    # P2P trading quantities: P2P injection and load with losses
    if options["grid"]:
        power_P2P_Loss = {}
        power_P2P_Loss["inj"] = {}
        power_P2P_Loss["load"] = {}
        power_P2P_Loss["inj"][z] = {}
        power_P2P_Loss["load"][z] = {}
        for j in nodes:
            power_P2P_Loss["inj"][z][j] = {}
            power_P2P_Loss["load"][z][j] = {}
            for t in time_steps:
                power_P2P_Loss["inj"][z][j][t] = model.addVar(vtype="C", name="power_from_P2P" + str(z) + "_t" + str(t))
                power_P2P_Loss["load"][z][j][t] = model.addVar(vtype="C", name="power_to_P2P" + str(z) + "_t" + str(t))

    # P2P trading quantities
    power_P2P_sum = {}
    power_P2P_sum["inj"] = {}
    power_P2P_sum["load"] = {}
    power_P2P_sum["inj"][z] = {}
    power_P2P_sum["load"][z] = {}

    for t in time_steps:
            power_P2P_sum["inj"][z][t] = model.addVar(vtype="C", name="power_from_P2P_sum" + str(z) + "_t" + str(t))
            power_P2P_sum["load"][z][t] = model.addVar(vtype="C", name="power_to_P2P_sum" + str(z) + "_t" + str(t))


    if options["grid"]:
        power_P2P_sum["Injloss"] = {}
        power_P2P_sum["Loadloss"] = {}
        power_P2P_sum["Injloss"][z] = {}
        power_P2P_sum["Loadloss"][z] = {}
        for t in time_steps:
            power_P2P_sum["Injloss"][z][t] = model.addVar(vtype="C", name="power_to_P2P_sum" + str(z) + "_t" + str(t))
            power_P2P_sum["Loadloss"][z][t] = model.addVar(vtype="C", name="power_to_P2P_sum" + str(z) + "_t" + str(t))

    # Residual network demand
    residual = {}
    residual["power"] = {}  # Residual network electricity demand
    residual["feed"] = {}  # Residual feed in
    for t in time_steps:
        residual["power"][t] = model.addVar(vtype="C", name="residual_power_t" + str(t))
        residual["feed"][t] = model.addVar(vtype="C", name="residual_feed_t" + str(t))



    # %% BALANCING UNIT VARIABLES

    # Electrical power to/from devices
    power = {}
    for device in ["from_grid", "to_grid", "PV"]:
        power[device] = {}
        for t in time_steps:
            power[device][t] = model.addVar(vtype="C", lb=0, name="power_" + device + "_t" + str(t))


    # Total operational costs
    operational_costs = model.addVar(vtype="C", lb=-gp.GRB.INFINITY, name="total_operational_costs")

    # Objective function
    obj = model.addVar(vtype="C", lb=-gp.GRB.INFINITY, name="obj")

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # DEFINE OBJECTIVE FUNCTION
    model.update()
    model.setObjective(obj)
    model.ModelSense = gp.GRB.MINIMIZE

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # ADD CONSTRAINTS

    # Device generation <= device capacity

    for t in time_steps:
            for device in ["HP", "EB"]:
                model.addConstr(heat_dom[device][z][t] <= nodes[z][device]["cap"], name=str(device) + "_cap_" + str(z))


    for t in time_steps:
            # Electric boiler
            model.addConstr(heat_dom["EB"][z][t] == power_dom["EB"][z][t],
                            name="EB_energybalance_" + str(z) + "_" + str(t))
            # Heat Pump
            model.addConstr(heat_dom["HP"][z][t] == power_dom["HP"][z][t] * nodes[z]["COP_sh"][:, idx][t],
                            name="heatpump_energybalance_heating" + str(z) + "_" + str(t))

    # constraint photolvotaics
    for t in time_steps:
            model.addConstr(power_dom["PV"][z][t] <= nodes[z]["pv_power"][:, idx][t],
                            name="max_pv_power_" + str(z) + "_" + str(t))

    # SOC of storages <= storage capacity

    for device in ["TES"]:
            for t in time_steps:
                model.addConstr(soc_dom[device][z][t] <= soc_nom["TES"][z],
                                name="max_cap_" + str(device) + "_" + str(z) + "_" + str(t))
                model.addConstr(soc_dom[device][z][t] >= nodes[z]["TES"]["min_soc"] * soc_nom["TES"][z],
                                name="min_cap_" + str(device) + "_" + str(z) + "_" + str(t))

    # battery constraints
    for t in time_steps:
            model.addConstr(ch_dom["BAT"][z][t] <= nodes[z]["BAT"]["max_ch"],
                            name="max_ch_cap_bat_" + str(z) + "_" + str(t))
            model.addConstr(dch_dom["BAT"][z][t] <= nodes[z]["BAT"]["max_dch"],
                            name="max_dch_cap_bat_" + str(z) + "_" + str(t))



            model.addConstr(soc_dom["BAT"][z][t] <= nodes[z]["BAT"]["max_soc"] * nodes[z]["BAT"]["cap"],
                            name="max_soc_bat_" + str(z) + "_" + str(t))

            model.addConstr(soc_dom["BAT"][z][t] >= nodes[z]["BAT"]["min_soc"] * nodes[z]["BAT"]["cap"],
                            name="max_soc_bat_" + str(z) + "_" + str(t))

    # P2P Constraints
    for j in nodes:
            for t in time_steps:
                if z == j:
                    model.addConstr(power_P2P["inj"][z][j][t] == 0,
                                    name="P2P_Constraint" + str(z) + "_" + str(j) + "_" + str(t))
                    model.addConstr(power_P2P["load"][z][j][t] == 0,
                                    name="P2P_Constraint" + str(z) + "_" + str(j) + "_" + str(t))


    # sum up injection and load at P2P
    for t in time_steps:
        model.addConstr(power_P2P_sum["inj"][z][t] == sum(power_P2P["inj"][z][j][t] for j in nodes), name="P2P" + str(z)  + "_" + str(t))
        model.addConstr(power_P2P_sum["load"][z][t] == sum(power_P2P["load"][z][j][t] for j in nodes), name="P2P" + str(z)  + "_" + str(t))

        model.addConstr(power_P2P_sum["Injloss"][z][t] == sum(power_P2P_Loss["inj"][z][j][t] for j in nodes),name="P2P" + str(z) + "_" + str(t))
        model.addConstr(power_P2P_sum["Loadloss"][z][t] == sum(power_P2P_Loss["load"][z][j][t] for j in nodes),name="P2P" + str(z) + "_" + str(t))

    if options["grid"]:
    # losses at P2P feed and power
        if options["grid"]:
            for j in nodes:
                for t in time_steps:
                    if Ubergabedaten["P_P2P_inj"][z][j][t] > 0:
                        model.addConstr(power_P2P_Loss["inj"][z][j][t] == power_P2P["inj"][z][j][t]*(1 - (Ubergabedaten["Verl"][z][j][t]/Ubergabedaten["P_P2P_inj"][z][j][t])))
                    else:
                        model.addConstr(power_P2P_Loss["inj"][z][j][t] == power_P2P["inj"][z][j][t])
            for j in nodes:
                for t in time_steps:
                    if Ubergabedaten["P_P2P_load"][z][j][t] > 0:
                        model.addConstr(power_P2P_Loss["load"][z][j][t] == power_P2P["load"][z][j][t]*(1 - (Ubergabedaten["Verl"][j][z][t] / Ubergabedaten["P_P2P_inj"][j][z][t])))
                    else:
                        model.addConstr(power_P2P_Loss["load"][z][j][t] == power_P2P["load"][z][j][t] )


    # secondary condition to avoid simultaneously feed and power
    if counter >= 3:
        for t in time_steps:
            if Netzdaten["P_inj_Ges"][z][t] >  Netzdaten["P_load_Ges"][z][t]:
               model.addConstr(res_dom["power"][z][t] == 0)
            elif Netzdaten["P_inj_Ges"][z][t] <  Netzdaten["P_load_Ges"][z][t]:
                model.addConstr(res_dom["feed"][z][t] == 0)


    # split up total consumption and feed into grid and P2P part
    for t in time_steps:
        model.addConstr(res_dom["power"][z][t] == power_P2P_sum["load"][z][t] + power_grid["load"][z][t])
        model.addConstr(res_dom["feed"][z][t] == power_P2P_sum["inj"][z][t] + power_grid["inj"][z][t])

        if options["grid"]:
            model.addConstr(res_dom["InjwithLoss"][z][t] == power_P2P_sum["Injloss"][z][t] + power_grid["inj"][z][t])
            model.addConstr(res_dom["LoadwithLoss"][z][t] == power_P2P_sum["Loadloss"][z][t] + power_grid["load"][z][t])




    # EV CONSTRAINTS CONSTRAINTS
    device = "EV"
    if options["ev_status"]:

            # Energy balance
            for t in time_steps:
                if t == 0:
                    soc_prev = soc_init[device][z]

                else:
                    soc_prev = soc_dom[device][z][t - 1]

                model.addConstr(soc_dom[device][z][t] == soc_prev
                                + ch_dom[device][z][t] * nodes[z][device]["eta_ch_ev"] * dt
                                - dch_dom[device][z][t] / nodes[z][device]["eta_dch_ev"] * dt -
                                nodes[z]["ev_dem_leave"][t])

                if t == last_time_step:
                    model.addConstr(soc_dom[device][z][t] == soc_init[device][z])

                model.addConstr(dch_dom[device][z][t] <= nodes[z]["ev_avail"][t] * nodes[z][device]["max_dch_ev"])
                model.addConstr(ch_dom[device][z][t] <= nodes[z]["ev_avail"][t] * nodes[z][device]["max_ch_ev"])


                model.addConstr(soc_dom[device][z][t] >= nodes[z][device]["min_soc"] * nodes[z][device]["cap"])
                model.addConstr(soc_dom[device][z][t] <= nodes[z][device]["max_soc"] * nodes[z][device]["cap"])

    else:

            for t in time_steps:
                model.addConstr(soc_dom[device][z][t] == 0)
                model.addConstr(dch_dom[device][z][t] == 0)
                model.addConstr(ch_dom[device][z][t] == 0)

    # %% DOMESTIC FLEXIBILITIES

    # SOC coupled over all times steps (Energy amount balance, kWh)
    # Energy balance TES
    device = "TES"
    for t in time_steps:
            if t == 0:
                soc_prev = soc_init[device][z]
            else:
                soc_prev = soc_dom[device][z][t - 1]

            model.addConstr(soc_dom[device][z][t] == soc_prev * nodes[z][device]["eta_tes"]
                            + ch_dom[device][z][t] * dt
                            - dch_dom[device][z][t] * dt, name="TES_storage_balance_" + str(z) + "_" + str(t))

            model.addConstr(ch_dom[device][z][t] == heat_dom["HP"][z][t], name="Heat_charging_" + str(z) + "_" + str(t))
            model.addConstr(dch_dom[device][z][t] <= nodes[z]["heat"][:, idx][t] + 0.5 * nodes[z]["dhw"][:, idx][t],
                            name="Heat_discharging_" + str(z) + "_" + str(t))

            if t == last_time_step:
                model.addConstr(soc_dom[device][z][t] == soc_init[device][z],
                                name="End_TES_Storage_" + str(z) + "_" + str(t))

    # Energy balance battery
    device = "BAT"
    for t in time_steps:
            if t == 0:
                soc_prev = soc_init[device][z]
            else:
                soc_prev = soc_dom[device][z][t - 1]

            model.addConstr(soc_dom[device][z][t] == soc_prev
                            + ch_dom[device][z][t] * nodes[z][device]["eta_bat"] * dt
                            - dch_dom[device][z][t] / nodes[z][device]["eta_bat"] * dt,
                            name="BAT_storage_balance_" + str(z) + "_" + str(t))

            if t == last_time_step:
                model.addConstr(soc_dom[device][z][t] == soc_init[device][z],
                                name="End_BAT_Storage_" + str(z) + "_" + str(t))

    # ENERGY BALANCES (Power balance, kW)

    for t in time_steps:
        # energy balance with losses
        if options["grid"]:

            model.addConstr(res_dom["LoadwithLoss"][z][t] + power_dom["PV"][z][t] + dch_dom["BAT"][z][t] + dch_dom["EV"][z][t]
                            == nodes[z]["elec"][:, idx][t] + power_dom["HP"][z][t] + power_dom["EB"][z][t] +
                            ch_dom["BAT"][z][t] + ch_dom["EV"][z][t] + res_dom["feed"][z][t],
                            name="Electricity_balance_" + str(z) + "_" + str(t))
        # energy balance without losses
        else:
            model.addConstr(
                res_dom["power"][z][t] + power_dom["PV"][z][t] + dch_dom["BAT"][z][t] + dch_dom["EV"][z][t]
                == nodes[z]["elec"][:, idx][t] + power_dom["HP"][z][t] + power_dom["EB"][z][t] +
                ch_dom["BAT"][z][t] + ch_dom["EV"][z][t] + res_dom["feed"][z][t],
                name="Electricity_balance_" + str(z) + "_" + str(t))

        model.addConstr(res_dom["feed"][z][t] <= power_dom["PV"][z][t] + dch_dom["BAT"][z][t] + dch_dom["EV"][z][t],
                            name="Feed-in_max_" + str(z) + "_" + str(t))


    for t in time_steps:
            # Heat balance
            model.addConstr(heat_dom["EB"][z][t] == 0.5 * nodes[z]["dhw"][:, idx][t],
                            name="Heat_balance_" + str(z) + "_" + str(t))


    # OBJECTIVE FUNCTIONS

    if options["grid"]:
        # Users operational costs with losses
        model.addConstr(operational_costs >= sum(power_grid["load"][z][t] * params["eco"]["elec_price"] - power_grid["inj"][z][t] * params["eco"]["eeg_pv"] -(power_P2P_sum["Injloss"][z][t])* admm_params["lambda_load"][z][t] + sum((power_P2P_Loss["load"][z][j][t]* admm_params["lambda_load"][j][t]) + (rho[counter] / 2) * (((power_P2P["inj"][z][j][t] - Ubergabedaten["P_P2P_load"][j][z][t]) ** 2) + ((Ubergabedaten["P_P2P_inj"][j][z][t] - power_P2P["load"][z][j][t]) ** 2))for j in nodes) for t in time_steps))

    else:


        # Users operational costs
        model.addConstr(
            operational_costs >= sum(
                power_grid["load"][z][t] * params["eco"]["elec_price"] - power_grid["inj"][z][t] * params["eco"][
                    "eeg_pv"] - power_P2P_sum["inj"][z][t] * admm_params["lambda_load"][z][t] + sum(
                    (power_P2P["load"][z][j][t] * admm_params["lambda_load"][j][t]) + (rho[counter] / 2) * (
                                ((power_P2P["inj"][z][j][t] - Ubergabedaten["P_P2P_load"][j][z][t]) ** 2) + (
                                    (Ubergabedaten["P_P2P_inj"][j][z][t] - power_P2P["load"][z][j][t]) ** 2))
                    for j in nodes) for t in time_steps))






    model.addConstr(obj == operational_costs)

    # START OPTIMIZATION
    # set objective function
    model.setObjective(obj, gp.GRB.MINIMIZE)           
    
    # adjust gurobi settings
    #model.Params.TimeLimit = 3600
    model.Params.NonConvex = 2
    model.Params.NumericFocus = 3
    model.Params.MIPGap = 0.05
    
    model.optimize()

    optiparams = {"runtime": model.runtime,
                  #"mipgap": model.MIPGap
                 }
    
    if model.status==gp.GRB.Status.INFEASIBLE:
        model.computeIIS()
        f=open('errorfile_stoch_central.txt','w')
        f.write('\nThe following constraint(s) cannot be satisfied:\n')
        for c in model.getConstrs():
            if c.IISConstr:
                f.write('%s' % c.constrName)
                f.write('\n')
        f.close()
    
    model.printStats()
    #%% RETRIEVE RESULTS

    if model.solCount == 0:
        print("Model is infeasible")
        model.computeIIS()
        model.write("model_iis.ilp")

    later = time.time()
    difference = later - now
    print("********************************************")
    print("Model run time was " + str(difference) + " seconds")
    print("********************************************")


    costs = obj.X

    res = dict(p_load_grid=np.array([power_grid["load"][z][t].X for t in time_steps]),
               p_inj_grid=np.array([power_grid["inj"][z][t].X for t in time_steps]),
               p_load_P2P_Ges=np.array([power_P2P_sum["load"][z][t].X  for t in time_steps]),
               p_inj_P2P_Ges=np.array([power_P2P_sum["inj"][z][t].X for t in time_steps]),
                p_inj_P2P={j: np.array([power_P2P["inj"][z][j][t].X for t in time_steps])for j in nodes },
                p_load_P2P={j: np.array([power_P2P["load"][z][j][t].X for t in time_steps])for j in nodes },
                power={dev: np.array([power_dom[dev][z][t].X for t in time_steps]) for dev in ["HP", "EB", "PV"]},
                heat={dev: np.array([heat_dom[dev][z][t].X for t in time_steps]) for dev in ["HP", "EB"]},
                soc={dev: np.array([soc_dom[dev][z][t].X for t in time_steps]) for dev in ["TES", "BAT", "EV"]},
                ch={dev: np.array([ch_dom[device][z][t].X for t in time_steps]) for dev in ["TES", "BAT", "EV"]},
                dch={dev: np.array([dch_dom[device][z][t].X for t in time_steps]) for dev in ["TES", "BAT", "EV"]},
                p_use= np.array([res_dom["power"][z][t].X for t in time_steps]),
                p_sell= np.array([res_dom["feed"][z][t].X for t in time_steps]),
                PV_P = np.array([power_dom["PV"][z][t].X for t in time_steps]))


    return costs, res, optiparams, difference

