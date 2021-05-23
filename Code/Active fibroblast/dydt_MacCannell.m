function output=dydt_MacCannell(t,X,flag_ode)
%Main code from which to run simulations, either from the original model
%of MacChannell et al. (2007)

%%  Define model parameters
%physical constants
R=8314.5;  % millijoule_per_mole_kelvin (in membrane) - Ideal gas constant
T=295;     % kelvin (in membrane) - Absolute temperature
F=96487;   % coulomb_per_mole (in membrane) - Faraday constant
Rtonf=(R*T)/F;
C_fm=6.3;  % picoF (in membrane) - Membrane capacitance
% extracellular ionic concentrations
fK_o=5.4;  % millimolar (in standard_ionic_concentrations)
fNa_o=140.0;  % millimolar (in standard_ionic_concentrations)
% initital intracellular ionic concentrations
fK_i=129.4349; % millimolar (in standard_ionic_concentrations)
fNa_i=10;  % millimolar (in standard_ionic_concentrations)

%%   Constants of fibroblast Parameter
fE_K=-87;
fg_Kv=0.25;
fr_Kv=0;
fs_Kv=1;
fg_K1=0.482;
vol_i=0.005884;% 
fV_rev=-150.0;
fB=-200;
fK_mk=1.0;%mM3.0
fK_mNa=11.0;%mM
fG_bNa=0.0095;
fI_NaKoo=2.002;

%% Global Variable for Time
global tStep tArray
if t > tArray(tStep) % Roughly eliminates data from rejected time steps
    tStep = tStep + 1;
end
tArray(tStep) = t;

% give names to the state vector values
v_f = X(1);
fr_Kv = X(2);
s_Kv = X(3);

%  reversal potentials
fE_Na = Rtonf*log(fNa_o/fNa_i);
fE_K=Rtonf*log(fK_o/fK_i);

%  calculate Time- and Voltage- Dependent K+ Current  (pA/pF)
ft_r=20.3+138*exp(-((v_f+20)/25.9)^2);
ft_s=1574+5268*exp(-((v_f+23.0)/22.7)^2);
%Active 1 Model
fr_ioo = 1.0/(1.0+exp(-(v_f+20.0)/11));
fs_ioo = 1.0/(1.0+exp((v_f+23.0)/7));
dfr_Kv = (fr_ioo-fr_Kv)/ft_r;
ds_Kv = (fs_ioo-fs_Kv)/ft_s;
fI_Kv=fg_Kv*fr_Kv*fs_Kv*(v_f-fE_K);
global fI_Kv_store
fI_Kv_store(tStep) = fI_Kv;

%   Timei-Independent Inward-Rectifer K+ Current, I_K1(pA/pF)
falpha_K1=0.1/(1.0+exp(0.06*(v_f-fE_K-200)));
fbate_K1=(3.0*exp(0.0002*(v_f-fE_K+100))+exp(0.1*(v_f-fE_K-10)))/...
    (1.0+exp(-0.5*(v_f-fE_K)));
ref_fIK1=falpha_K1/(falpha_K1+fbate_K1);
fI_K1=fg_K1*ref_fIK1*(v_f-fE_K);
global fI_K1_store
fI_K1_store(tStep) = fI_K1;
%    Na+ - K+ Exchanger Current, I_NaK(pA/pF)
fI_NaK=fI_NaKoo*(fK_o/(fK_o+fK_mk))*((v_f-fV_rev)/(v_f-fB))*...
    ((fNa_i/(fNa_i+fK_mNa))^1.5);
global  fI_NaK_store
fI_NaK_store(tStep) = fI_NaK;
%Background Na+ current, Ib,Na(pA/pF)
fI_bNa=fG_bNa*(v_f-fE_Na);%
global fI_bNa_store
fI_bNa_store(tStep) = fI_bNa;

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

fI_tot = fI_Kv+fI_K1+fI_NaK+fI_bNa;
global fI_tot_store
fI_tot_store(tStep) = fI_tot;

dv_f = -fI_tot-Istim;

%output the state venctor when ode_flag==1, and the calculated currents and fluxes otherwise
if flag_ode==1
    output=[dv_f dfr_Kv ds_Kv]';
else
    output=[fI_Kv fI_K1 fI_NaK fI_bNa];
end

