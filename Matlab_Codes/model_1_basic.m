% AMT CH3030 Term Project
% Deepanjhan Das | CH22B020
% model 1 (simplified model)

% Full Simulation Code for Ethanol Dehydration by Pervaporation
% Mass and Energy Balances included
clear; clc;
close all;

%% Parameters
% Membrane and Process Parameters
membrane_thickness = 100e-6; % Membrane thickness in meters (100 microns)
area = 1.0; % Membrane area in square meters
P_water = 1e-9; % Permeability of water through PDMS (m^2/s)
P_ethanol = 1e-11; % Permeability of ethanol through PDMS (m^2/s)
R = 8.314; % Gas constant (J/mol·K)
Mw_water = 18.015; % Molecular weight of water (g/mol)
Mw_ethanol = 46.07; % Molecular weight of ethanol (g/mol)

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
    % Concentration and Temperature-dependent Permeabilities (Arrhenius equation)
    P_water_eff = P_water * exp(-500 / temperature_feed(i)); % Water permeability (arbitrary activation energy)
    P_ethanol_eff = P_ethanol * exp(-700 / temperature_feed(i)); % Ethanol permeability (arbitrary activation energy)

    % Calculate flux for water and ethanol (Fick's law)
    flux_water(i) = P_water_eff * (water_concentration_feed(i)) / membrane_thickness;
    flux_ethanol(i) = P_ethanol_eff * (ethanol_concentration_feed(i)) / membrane_thickness;

    % Mass balance for feed concentration changes
    % Using approximate feed volume change due to permeation
    d_water = -flux_water(i) * area * dt / feed_volume;
    d_ethanol = -flux_ethanol(i) * area * dt / feed_volume;

    % Update concentrations in feed
    water_concentration_feed(i+1) = water_concentration_feed(i) + d_water;
    ethanol_concentration_feed(i+1) = ethanol_concentration_feed(i) + d_ethanol;
    
    % Energy balance (Assume a constant heat loss coefficient)
    heat_capacity_feed = 4.18 * 1000; % Approx heat capacity of ethanol-water mix (J/kg·K)
    heat_loss_coefficient = 10; % Arbitrary heat loss coefficient
    dT = - (heat_loss_coefficient * (temperature_feed(i) - initial_temperature) * dt) / ...
         (feed_volume * 1000 * heat_capacity_feed); % Cooling effect in °C/s
    
    % Update temperature
    temperature_feed(i+1) = temperature_feed(i) + dT;

    % Ensure concentration values stay within logical limits
    if water_concentration_feed(i+1) < 0
        water_concentration_feed(i+1) = 0;
    end
    if ethanol_concentration_feed(i+1) < 0
        ethanol_concentration_feed(i+1) = 0;
    end
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
grid on;
plot(time/60, flux_ethanol, 'r-', 'LineWidth', 1.5);
xlabel('Time (min)');
ylabel('Flux (mol/m^2·s)');
legend('Water Flux', 'Ethanol Flux', Location='best');
title('Water and Ethanol Flux over Time');