%Main code from which to run simulations, either from the original model
%of Sachse et al. (2007)
% Link to Article:
% https://link.springer.com/article/10.1007/s10439-007-9405-8

clc
clear all
%initial conditions for state variables
V_f = -60.1923588449342;
C_0Shkr = 0.911000000000000;
C_1Shkr = 0.0857000000000000;
C_2Shkr = 0.00302000000000000;
C_3Shkr = 4.74000000000000e-05;
C_4Shkr = 4.74000000000000e-05;
O_Shkr = 0;

%X0 is the vector for initial sconditions for state variables
X0=[V_f C_0Shkr C_1Shkr C_2Shkr C_3Shkr C_4Shkr O_Shkr]';

CL=1000;%pacing cycle length in ms
beats=1000;%number of beats in the simulation

global tStep tArray
global  fI_Kir_store fI_Shkr_store fI_b_store Istim_store fI_tot_store

tStep = 1;
tArray = zeros(1,1e7);
fI_Kir_store = zeros(1,1e7);
fI_Shkr_store = zeros(1,1e7);
fI_b_store = zeros(1,1e7);
Istim_store = zeros(1,1e7);
fI_tot_store = zeros(1,1e7);
%% Run simulation
tic
tspan = [0 100e2];
options = odeset('RelTol',1e-5,'MaxStep',1);
[time,X] = ode15s(@dydt_Sachse,tspan,X0,options,1);
toc

% %% Output variables
% 
% for n=[1:beats]
%     [time X]=ode15s(@model,[0 CL],X0,options,1);
%     X0=X(size(X,1),:);
%     n; %output beat number to the screen to monitor runtime progress
% end
%% Output variables
tArray = tArray(1:tStep);
fI_Kir = fI_Kir_store(1:tStep);
fI_Shkr = fI_Shkr_store(1:tStep);
fI_b = fI_b_store(1:tStep);
fI_tot = fI_tot_store(1:tStep);
Istim = Istim_store(1:tStep);
fI_K = fI_Kir+fI_Shkr;
%rename values in the state variables vector
V_f=X(:,1);
fr_Kv=X(:,2);
fs_Kv=X(:,3);
%calculate and name dependent variables for the final beat in the
%simulation (i.e. currents and fluxes)

%create plots showing results for the final paced beat
% figure
hold on
subplot(2,1,1),plot(time,V_f),title('V_f')
title('Membrane potential, Vf(mV)','Fontsize',18);
hold on
subplot(2,1,2),plot(tArray,Istim),title('fI_tot')
title('Membrane potential, Vf(mV)','Fontsize',18);