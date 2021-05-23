function output=dydt_Sachse(t,X,flag_ode)
%Main code from which to run simulations, either from the original model
%of Sachse et al. (2008)

%%  Define model parameters
%physical constants
R=8314.5;  % millijoule_per_mole_kelvin (in membrane) - Ideal gas constant
T=295;     % kelvin (in membrane) - Absolute temperature
F=96487;   % coulomb_per_mole (in membrane) - Faraday constant
Rtonf=(R*T)/F;
C_fm = 4.5;  % picoF (in membrane) - Membrane capacitance
% extracellular ionic concentrations
fK_o = 5.4;  % millimolar (in standard_ionic_concentrations)
fNa_o = 140.0;  % millimolar (in standard_ionic_concentrations)
% initital intracellular ionic concentrations
fK_i = 129.4349; % millimolar (in standard_ionic_concentrations)
fNa_i = 10;  % millimolar (in standard_ionic_concentrations)

%%   Constants of fibroblast Parameter for [K_o]=5mM
P_Shkr=5.4e-6;
kvo=30e-3;%
k_vo=2e-3;%
zv=1.28;%
z_v=-1.53;%
ko=77e-3;%
k_o=18e-3;%

G_Kir=1.02;%
a_Kir=0.94;
b_Kir=1.26;

fG_b=6.9e-3;
fE_b=0.0;
%% Global Variable for Time
global tStep tArray
if t > tArray(tStep) % Roughly eliminates data from rejected time steps
    tStep = tStep + 1;
end
tArray(tStep) = t;

% give names to the state vector values
V_f = X(1);
C_0Shkr = X(2);
C_1Shkr = X(3);
C_2Shkr = X(4);
C_3Shkr = X(5);
C_4Shkr = X(6);
O_Shkr = X(7);


%  reversal potentials
fE_Na = Rtonf*log(fNa_o/fNa_i);
fE_K = Rtonf*log(fK_o/fK_i);

%%  I_Kir
O_Kir=1/(a_Kir+exp(b_Kir*(V_f-fE_K)/Rtonf));
fI_Kir=G_Kir*O_Kir*sqrt(fK_o*0.001)*(V_f-fE_K);
global fI_Kir_store
fI_Kir_store(tStep) = fI_Kir;

%%  I_Shkr
kv=kvo*exp(V_f*zv/Rtonf);
k_v=k_vo*exp(V_f*z_v/Rtonf);
if C_0Shkr<0
    C_0Shkr=0;
    dC_0Shkr = 0;
elseif C_0Shkr>1
    C_0Shkr=1;
    dC_0Shkr = 1;
else
    dC_0Shkr=(k_v*C_1Shkr-4*kv*C_0Shkr);
end
if C_1Shkr<0
    C_1Shkr=0;
    dC_1Shkr =0;
elseif C_1Shkr>1
    C_1Shkr=1;
    dC_1Shkr = 1;
else
    dC_1Shkr=(2*k_v*C_2Shkr+4*kv*C_0Shkr-(3*kv+k_v)*C_1Shkr);
end
if C_2Shkr<0
    C_2Shkr=0;
    dC_2Shkr = 0;
elseif C_2Shkr>1
    C_2Shkr=1;
    dC_2Shkr = 0;
else
    dC_2Shkr=(3*k_v*C_3Shkr+3*kv*C_1Shkr-(2*kv+2*k_v)*C_2Shkr);
end
if C_3Shkr<=0
    C_3Shkr=0;
    dC_3Shkr = 0;
elseif C_3Shkr>=1
    C_3Shkr=1;
    dC_3Shkr = 0;
else
    dC_3Shkr=(4*k_v*C_4Shkr+2*kv*C_2Shkr-(kv+3*k_v)*C_3Shkr);
end
if C_4Shkr<=0
    C_4Shkr=0;
    dC_4Shkr = 0;
elseif C_4Shkr>=1
    C_4Shkr=1;
    dC_4Shkr = 0;
else
    dC_4Shkr=(k_o*O_Shkr+kv*C_3Shkr-(ko+4*k_v)*C_4Shkr);
end
if O_Shkr<0
    O_Shkr=0;
    dO_Shkr=0;
elseif O_Shkr>1
    O_Shkr=1;
    dO_Shkr =0;
else
    dO_Shkr=(ko*C_4Shkr-k_o*O_Shkr);
end
fI_Shkr=P_Shkr*O_Shkr*(V_f*F/Rtonf)*(fK_i-fK_o*exp(-V_f/Rtonf))/...
    (1-exp(-V_f/Rtonf));
global fI_Shkr_store
fI_Shkr_store(tStep) = fI_Shkr;

%%  I_b & I_stim
fI_b=fG_b*(V_f-fE_b);
global fI_b_store
fI_b_store(tStep) = fI_b;

%  simulation time
cycleLength = 1000;
stimstrength = -20;
if mod(t-20,cycleLength) <= 1
    Istim = stimstrength;
else
    Istim = 0.0;
end
global Istim_store
Istim_store(tStep) = Istim;

fI_tot = fI_Kir + fI_Shkr + fI_b;
global fI_tot_store
fI_tot_store(tStep) = fI_tot;

dV_f = -fI_tot-Istim;

%output the state venctor when ode_flag==1, and the calculated currents and fluxes otherwise
if flag_ode==1
    output=[dV_f dC_0Shkr dC_1Shkr dC_2Shkr dC_3Shkr dC_4Shkr dO_Shkr]';
else
    output=[fI_Kv fI_K1 fI_NaK fI_bNa];
end

