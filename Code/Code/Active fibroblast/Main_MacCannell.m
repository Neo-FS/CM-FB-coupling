%Main code from which to run simulations, either from the original model
%of MacChannell et al. (2007)
% Link to Article:
% https://www.sciencedirect.com/science/article/pii/S0006349507712092

clc
clear all
%initial conditions for state variables
v_f= -47.9595056382552;
fr_Kv = 0.0729811147816670;
fs_Kv= 0.906050007963956;
%X0 is the vector for initial sconditions for state variables
X0=[v_f fr_Kv fs_Kv]';

CL=1000;%pacing cycle length in ms
beats=1000;%number of beats in the simulation

global tStep tArray
global  fI_Kv_store fI_K1_store fI_NaK_store fI_bNa_store Istim_store fI_tot_store

tStep = 1;
tArray = zeros(1,1e7);
fI_Kv_store = zeros(1,1e7);
fI_K1_store = zeros(1,1e7);
fI_NaK_store = zeros(1,1e7);
fI_bNa_store = zeros(1,1e7);
Istim_store = zeros(1,1e7);
fI_tot_store = zeros(1,1e7);
%% Run simulation
tic
tspan = [0 100e2];
options = odeset('RelTol',1e-5,'MaxStep',1);
[time,X] = ode15s(@dydt_MacCannell,tspan,X0,options,1);
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
fI_Kv = fI_Kv_store(1:tStep);
fI_K1 = fI_K1_store(1:tStep);
fI_NaK = fI_NaK_store(1:tStep);
fI_bNa = fI_bNa_store(1:tStep);
Istim = Istim_store(1:tStep);
fI_tot = fI_tot_store(1:tStep);
fI_K = fI_Kv+fI_K1;
%rename values in the state variables vector
v_f=X(:,1);
fr_Kv=X(:,2);
fs_Kv=X(:,3);
%calculate and name dependent variables for the final beat in the
%simulation (i.e. currents and fluxes)

%create plots showing results for the final paced beat
% figure
hold on
subplot(2,1,1),plot(time,v_f),title('V_f')
title('Membrane potential, Vf(mV)','Fontsize',18);
hold on
subplot(2,1,2),plot(tArray,Istim),title('fI_tot')
title('Stimulation current, Istim (mV)','Fontsize',18);