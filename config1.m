function [x,b,cons] = config1()

secyr = 3600*24*365; % seconds per year

% Material parameters
cons.A = 1e-24*secyr; % Glen's Flow law constant [Pa^{-3}*yr^{-1}]
cons.g = 9.81; % gravitational constant [m/s^2]
cons.rho = 910.0; % density of ice [kg/m^3]

% Debris characteristics
cons.D0 = 0.05; % characteristic debris thickness [m]
cons.c = 0.000; % debris concentration

% Boundary condition
cons.H_t = 30;% terminal ice cliff height [m]

% Surface mass balance forcing
cons.amax = 2; % surface mass balance (SMB) cutoff [m/yr]
cons.ELA = 3000; % equilibrium line altitute [m]
cons.gamma = 0.0075; % SMB gradient 

% Terminal cryokarst
cons.cryo = 0; % boolean - terminal cryokart
cons.lambdam = 0.2; % maximum cryokarst area fraction
cons.maxth = 1.1e5; % cryokarst upper driving stress threshold [Pa] 
cons.minth = 0.6e5; % cryokarst lower driving stress threshold [Pa]

% Computational grid

L = 14000; % Length [m]
dx = 25; % grid spacing [m]
x = [0:dx:L]';

% Bed - should be the same size as the grid

% This bed has a short steep headwall followed by a less steep main bed
theta = 0.1; 
b1 = 4000-x;
b2 = 3000-theta*(x-1000);
b = b1;
b(x>1000) = b2(x>1000);








