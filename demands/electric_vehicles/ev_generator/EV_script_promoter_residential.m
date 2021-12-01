close all;
clearvars;
clc;
pattern = '_promoter_commercial';
weekend = '';
ev_data = ['EV_data' pattern weekend '.mat'];
load(ev_data);

%% input parameters
num_profiles = 550;     % number of profiles
charge_power = 43.0;      % chargin power in kW
year = 2017;            % year of simulation
charger_eff = 0.97;     % charger efficiency

% from EV_data.mat:
% type_car: e-vehilce type acc. to Fraunhofer study, see type_description for
% parameters; 
% small: small car; medium: medium car; big: big car; other: transporter

type ='small';

if strcmp(type,'small')
    bat_cap = type_car(1,1);        % battery capacity in kWh      
    spec_demand = type_car(2,1);    % spec. demand in kWh/km
elseif strcmp(type,'medium')
    bat_cap = type_car(1,2);        % battery capacity in kWh      
    spec_demand = type_car(2,2);    % spec. demand in kWh/km
elseif strcmp(type,'big')
    bat_cap = type_car(1,3);        % battery capacity in kWh      
    spec_demand = type_car(2,3);    % spec. demand in kWh/km
elseif strcmp(type,'other')
    bat_cap = type_car(1,3);        % battery capacity in kWh      
    spec_demand = type_car(2,3);    % spec. demand in kWh/km
end
    
%% stochastic parameters

% from EV_data.mat:
% Arrival_prob: from Mobilität in Deutschland 2008 study
% Departure_prob: from Mobilität in Deutschland 2008 study
% Daily_distance_prob: from Diss. Michael Agsten
% Distance_prob_factor: from Mobilität in Deutschland 2008 study
% Velocity_prob: from Diss. Lin Zhao

%% open parallel pool (Matlab R2013)
% if matlabpool('size')==0
%     matlabpool open    
% end

%% open parallel pool (Matlab R2015)
p = gcp('nocreate'); 

%% run stoch. EV-profiles 

[EV_occurrence,EV_daily_demand,EV_uncontrolled_charging_profile, ...
    SOC_uncontrolled_charging_profile]=EV_generator(num_profiles, ...
charge_power, year, charger_eff, bat_cap, spec_demand, Arrival_prob, ...
Departure_prob,Daily_distance_prob,Distance_prob_factor,Velocity_prob); 

%% output

% EV_occurrence: occurrence profile for EV at home (1: at home; 0 not at
% home)
% EV_daily_demand: daily charging demand in kWh
% EV_uncontrolled_charging_profile: power profile for uncontrolled charging
% in kW
% SOC_uncontrolled_charging_profile: state of charge (SOC) in kWh

if ~exist ('ev_profiles','dir')
    mkdir('ev_profiles')
end

save([pwd '\ev_profiles\' 'EV' pattern '_' type '_' num2str(charge_power) '_' num2str(year) '.mat'],'EV_occurrence','EV_daily_demand', ...
'EV_uncontrolled_charging_profile','SOC_uncontrolled_charging_profile')

%% close parallel pool (Matlab R2013)
%matlabpool close

%% close parallel pool (Matlab R2015)
%delete(p);
