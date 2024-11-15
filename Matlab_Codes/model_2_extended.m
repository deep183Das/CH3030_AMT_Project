% AMT CH3030 Term Project
% Deepanjhan Das | CH22B020
% model 2 (extension of model 1)

% Full Simulation Code for Ethanol Dehydration by Pervaporation
% Extended Mass and Energy Balances with VLE and Membrane Fouling
clear; clc;
close all;

%% Parameters
% Membrane and Process Parameters
membrane_thickness = 100e-6; % Membrane thickness in meters (100 microns)
area = 1.0; % Membrane area in square meters
P_water_0 = 1e-9; % Initial permeability of water through PDMS (m^2/s)
P_ethanol_0 = 1e-11; % Initial permeability of ethanol through PDMS (m^2/s)
R = 8.314; % Gas constant (J/mol·K)
Mw_water = 18.015; % Molecular weight of water (g/mol)
Mw_ethanol = 46.07; % Molecular weight of ethanol (g/mol)
H_vap_water = 40650; % Heat of vaporization for water (J/mol)
H_vap_ethanol = 38100; % Heat of vaporization for ethanol (J/mol)
fouling_factor = 0.999; % Membrane fouling factor per time step

% Feed conditions
initial_concentration_water = 0.1; % Initial mole fraction of water in feed
initial_temperature = 343; % Temperature in Kelvin (70 °C)
feed_volume = 0.5; % Volume of feed in liters

% Simulation parameters
time_end = 3600; % Simulation time in seconds
dt = 10; % Time step in seconds

% Arrays to store results
time = 0:dt:time_end;
water_concentration_feed = zeros(size(time));
ethanol_concentration_feed = zeros(size(time));
temperature_feed = zeros(size(time));
flux_water = zeros(size(time));
flux_ethanol = zeros(size(time));

% Initial conditions
water_concentration_feed(1) = initial_concentration_water;
ethanol_concentration_feed(1) = 1 - initial_concentration_water;
temperature_feed(1) = initial_temperature;

%% Simulation Loop
for i = 1:length(time)-1
    % Update permeabilities with fouling effect
    P_water = P_water_0 * (fouling_factor^i);
    P_ethanol = P_ethanol_0 * (fouling_factor^i);

    % Concentration and Temperature-dependent Permeabilities (Arrhenius equation)
    P_water_eff = P_water * exp(-500 / temperature_feed(i)); % Water permeability
    P_ethanol_eff = P_ethanol * exp(-700 / temperature_feed(i)); % Ethanol permeability

    % VLE using Raoult's Law (Assuming ideal solution behavior for simplicity)
    P_sat_water = exp(16.3872 - 3885.7 / (temperature_feed(i) - 42.98)); % Water saturation vapor pressure (Pa)
    P_sat_ethanol = exp(16.8958 - 3795.17 / (temperature_feed(i) - 42.98)); % Ethanol saturation vapor pressure (Pa)

    partial_pressure_water = P_sat_water * water_concentration_feed(i);
    partial_pressure_ethanol = P_sat_ethanol * ethanol_concentration_feed(i);

    % Flux calculation with VLE effect
    flux_water(i) = P_water_eff * (partial_pressure_water / membrane_thickness);
    flux_ethanol(i) = P_ethanol_eff * (partial_pressure_ethanol / membrane_thickness);

    % Mass balance for feed concentration changes
    d_water = -flux_water(i) * area * dt / feed_volume;
    d_ethanol = -flux_ethanol(i) * area * dt / feed_volume;

    % Update concentrations in feed
    water_concentration_feed(i+1) = water_concentration_feed(i) + d_water;
    ethanol_concentration_feed(i+1) = ethanol_concentration_feed(i) + d_ethanol;
    
    % Energy balance considering heat of vaporization
    heat_capacity_feed = 4.18 * 1000; % Approx heat capacity of ethanol-water mix (J/kg·K)
    heat_loss_coefficient = 10; % Arbitrary heat loss coefficient
    dT = (- (flux_water(i) * H_vap_water + flux_ethanol(i) * H_vap_ethanol) * area * dt - ...
          heat_loss_coefficient * (temperature_feed(i) - initial_temperature) * dt) / ...
          (feed_volume * 1000 * heat_capacity_feed);

    % Update temperature
    temperature_feed(i+1) = temperature_feed(i) + dT;

    % Ensure concentration values stay within logical limits
    water_concentration_feed(i+1) = max(water_concentration_feed(i+1), 0);
    ethanol_concentration_feed(i+1) = max(ethanol_concentration_feed(i+1), 0);
end

%% Plot Results
figure;
subplot(3,1,1);
plot(time/60, water_concentration_feed, 'b-', 'LineWidth', 1.5);
xlabel('Time (min)');
ylabel('Water Concentration in Feed');
grid on;
title('Water Concentration in Feed vs Time');

subplot(3,1,2);
plot(time/60, ethanol_concentration_feed, 'r-', 'LineWidth', 1.5);
xlabel('Time (min)');
ylabel('Ethanol Concentration in Feed');
grid on;
title('Ethanol Concentration in Feed vs Time');

subplot(3,1,3);
plot(time/60, temperature_feed, 'k-', 'LineWidth', 1.5);
xlabel('Time (min)');
ylabel('Temperature in Feed (K)');
grid on;
title('Temperature in Feed vs Time');

figure;
plot(time/60, flux_water, 'b-', 'LineWidth', 1.5);
hold on;
plot(time/60, flux_ethanol, 'r-', 'LineWidth', 1.5);
xlabel('Time (min)');
ylabel('Flux (mol/m^2·s)');
legend('Water Flux', 'Ethanol Flux', Location='best');
grid on;
title('Water and Ethanol Flux over Time');