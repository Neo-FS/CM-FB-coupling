% function [M1]= Pandit_MacCannel
% MODEL NAME: Pandit_MacCannel in rat
% SHORT DESCRIPTION: This model reproduces the action potential recorded experimentally
% for myocytes isolated from the adult rat left ventricle (Pandit et al. 2001) and
% active fibroblasts from MacCannell (MacCannell et al. 2007).
% Geometrical structure of left ventricular myocytes in rat are redesign
% baed on the the experimental measeuremnt (Wulfsohn et al, The Anatomical
% Record Part A: Discoveries in Molecular, Cellular, and Evolutionary
% Biology, 2004) and Luo report (Luo CH and  Rudy Y. A, Circulation Research. 1994)
% Becasue Pandit model cannot simulate the Ca2+ cycling of rat myocytes effectively,
% we replaced the Pandit model formulation for Ca2+ cycling with that of the
% ten Tusscher model (ten Tusscher et al. 2006).Additionally, Maximal Irel
% conductance, V_rel, is adjusted to 0.102 mM/ms, a transition rate of I_rel, k4
% is adjust to 0.005 ms-1,to agree well with the experimental data (Kaprielian et al, 1999)
% Maximum INaK current,I_NaK_max, and Scaling factor for INaCa , K_NaCa, are adjusted to
% satisfy the convergencein of Na+，K+ and Ca2+ in epicardial myocytes and endocardial myocytes.

clc; clear all
t1=clock;
%------------------------------------------------------------------------
%                  Parameter for rat ventricular myocytes
%------------------------------------------------------------------------
Num_Fb = 0;  % the numeb of active fibroblast
Celltype = 1; % 1,Epicardial cell; 2,Endocardial cell
%  Physical constant
R = 8314.5;
T = 295;
F = 96487;
Rtonf = (R*T)/F;
% Intracellular volumes
Cm = 0.132;% microF

% Geometrical structure of cardiomyocytes
Length = 0.11; width = 0.028; depth = width/3; % millimeter
Cell_volume = 1000*pi*Length*width*depth/4;%microliter
Vc = Cell_volume*0.68;%9.36e-3; %microliter-Myoplasm volume
Vsr = 0.0667*Vc; % microliter-Sarcoplasmic reticulum volume
Vss = 0.0033*Vc; % microliter-Submembrance volume
inverseVcF = 1/(Vc*F);
inverseVcF2 = 1/(2*Vc*F);
inversevssF2 = 1/(2*Vss*F);

%External concentrations
Ca_o = 1.2;
Na_o=140;
K_o = 5.4;
% Initial internal concentrations
Na_i=10.737388188297878;
K_i = 1.392770305344122e+02;
Ca_i = 1.130719597532137e-04;
Ca_ss = 1.615093898432590e-04;
Ca_SR = 3.064016252685225;

% Calcium buffering dynamics
Bufc=0.20;
Kbufc=0.001;
Bufsr=10.;
Kbufsr=0.3;
Bufss=0.40;
Kbufss=2.5e-4;
Vmaxup=6.375e-3;
Kup = 2.5e-4;
Vrel = 0.102;%%40.8;%%
k1_= 0.15;
k2_= 0.045;
k3 = 0.060;
k4 = 0.005;%%1.5e-5;%%
EC=1.5;
maxsr = 2.5;
minsr = 1.;
Vleak = 3.6e-4;
Vxfer = 3.8e-3;

%   real Cm nanoF;
Cm=0.1;

%  Parameter of rat myocytes
if Celltype == 1
    G_Na=0.8;
    G_t = 0.035;   %   G_t microS;
    a = 0.886;  %   real a dimensionless;
    b = 0.114;
    K_NaCa=0.000009984*0.6482;
    I_NaK_max=0.08*0.962;
elseif Celltype == 2
    G_Na=0.8*1.33;
    G_t = 0.035*0.4547;   %   G_t microS;
    a = 0.583;  %   real a dimensionless;
    b = 0.417;
    K_NaCa=0.000009984*1.231;
    I_NaK_max=0.08*1.317;
else
    error('Error:Incorrect input of parameters')
end
G_Ca_L=0.031;
E_Ca_L=65;
tau_Ca_inact=9;
G_ss=0.007;
G_K1=0.024;
G_f=0.00145;
f_Na=0.2;
G_bNa=0.00008015;
G_bCa=0.0000324;
G_bK=0.000138;
K_m_K=1.5;
K_m_Na=10;
I_pCa_max=0.004;
d_NaCa=0.0001;
gamma_NaCa=0.5;

% Initial value of ion channge gate in myocytes
V=-80.3073;
sodium_m= 0.0044;
sodium_h=0.6666;
sodium_j= 0.6665;
r=0.0022;
s=0.9848;
s_slow=0.6511;
y=0.0036;
r_ss=0.0030;
s_ss=0.3166;
d=0.0000;
f_11=1;
f_12=1;
Ca_inact=0.9840;
sRR= 0.9931;

%--------------------------------------------------------------------------
%              Parameter for MacCannell== active fibroblast model
%--------------------------------------------------------------------------
C_fm=6.3;  % picroF (in membrane) - active firboblast  membrane capacitance
fK_o=K_o; % millimoles per liter
fNa_o=Na_o; % millimoles per liter
fK_i=129.4349; % millimoles per liter
fNa_i=10; % millimoles per liter

G_gap=3.0;  % nanoS
fG_Kv = 0.25; % nanoS per picoF;
fG_K1 = 0.4822; % nanoS per picoF;
fG_bNa=0.0095; % nanoS per picoF;

fE_K=-87;  % milliVolt
fV_rev=-150.0; % milliVolt
fB=-200; % milliVolt
fK_mk=1.0;% millimoles per liter
fK_mNa=11.0;% millimoles per liter
fI_NaK_inf=2.002; % picoA per picoF;

% Initial internal variable of active fibroblast
v_f=-48.747871371110890;
fr_Kv=0.068278257548253;
fs_Kv=0.975354703085357;

%   time of the simulation
PD =[];  PI =[];cycl =[];
%--------------------------------------------------------------------------
%              Parameter for simulation duration
%--------------------------------------------------------------------------
xns = 10;
dt=0.005;
sdur = 1;
I_stim = 5;
stimstrength = -6;  %  stimstrength
stims1 = 1000;
stims2 = 1000;
tbegin=20;
xnstims1 = stims1+tbegin+sdur;%+tbegin;     ???
xnstims2 =stims1*(xns-1)+tbegin+stims2+sdur;%tbegin;
endtime=xnstims2+stims2;  %duration of the simulation
nswitch = 0;
ncounts1=1;
napdd=0;
apdtime = zeros(xns,1);
time=0; step=0; % timestep(ms)
st=[];sv=[];svf =[];sa=[];sb=[];sCa=[];sd=[];se=[];sf=[];sg=[];APV = [];
while time<=endtime
    %--------------------------------------------------------------------------
    %            Simulation protocols
    %--------------------------------------------------------------------------
    if (time >=tbegin&&time<=tbegin+sdur)
        I_stim = stimstrength;
    elseif (time >=xnstims1-sdur&&time<=xnstims1)
        I_stim = stimstrength;
        nswitch = 1;
    else
        I_stim = 0;
    end
    if (time>xnstims1&& ncounts1 < xns-1 && nswitch==1)
        ncounts1=ncounts1+1;
        xnstims1=ncounts1*stims1+tbegin+sdur;
        nswitch=0;
    end
    if (time >=xnstims2-sdur&&time<=xnstims2)
        I_stim = stimstrength;
        skip=1;
        ncount=0;
    end
    
    %--------------------------------------------------------------------------
    %             Computer myocytes current
    %--------------------------------------------------------------------------
    %  Reverse voltage
    E_Na = Rtonf*log(Na_o/Na_i);
    E_K = Rtonf*log(K_o/K_i);
    E_Ca = 0.5*Rtonf*log(Ca_o/Ca_i);
    
    % Fast I_Na_current
    m_inf=(1/(1+exp((V+45)/(-6.5))));
    tau_m=(1.36/(((0.32*(V+47.13))/(1-exp(((-1)*0.1)*(V+(47.13 )))))+(0.08*exp(((-1)*V)/(11 )))));
    sodium_m=sodium_m+dt*((m_inf-sodium_m)/tau_m);
    
    h_inf=(1/(1+exp((V+(76.1 ))/(6.07 ))));
    j_inf=(1/(1+exp((V+(76.1 ))/(6.07 ))));
    if V>=-40
        tau_h=(0.4537 )*(1+exp(((-1)*(V+(10.66 )))/(11.1 )));
        tau_j=((11.63 )*(1+exp(((-1)*0.1)*(V+(32 )))))/exp(((-1)*2.535E-7)*V);
    else
        tau_h=(3.49 )/(((0.135*exp(((-1)*(V+(80 )))/(6.8 )))+(3.56*exp(0.079*V)))+(310000*exp(0.35*V)));
        tau_j=(3.49 )/((((V+(37.78 ))/(1+exp(0.311*(V+(79.23 )))))*...
            ((((-1)*127140)*exp(0.2444*V))-(3.474E-5*exp(((-1)*0.04391)*V))))...
            +((0.1212*exp(((-1)*0.01052)*V))/(1+exp(((-1)*0.1378)*(V+(40.14 ))))));
    end
    sodium_h=sodium_h+dt*((h_inf-sodium_h)/tau_h);
    sodium_j=sodium_j+dt*((j_inf-sodium_j)/tau_j);
    
    I_Na=((((G_Na*(sodium_m^3))*sodium_h)*sodium_j)*(V-E_Na));
    
    % background_currents
    I_bNa=(G_bNa*(V-E_Na));
    I_bCa=(G_bCa*(V-E_Ca));
    I_bK=(G_bK*(V-E_K));
    I_b=I_bNa+I_bCa+I_bK;
    
    %  Ca_independent_transient_outward_K_currents
    r_inf=(1/(1+exp((V+(10.6))/((-1)*(11.42 )))));
    tau_r=((1000 )/((45.16*exp(0.03577*(V+(50 ))))+(98.9*exp(((-1)*0.1)*(V+(38 ))))));
    r=r+dt*((r_inf-r)/tau_r);
    
    if Celltype == 1
        tau_s=((350*exp((-1)*(((V+(70 ))/(15 ))^2)))+35);
        tau_s_slow=(((3700)*exp((-1)*(((V+(70))/(30))^2)))+(35));
    elseif Celltype == 2
        tau_s = 550*exp(-(((V+70)/25)^2))+49;
        tau_s_slow = 3300*exp(-(((V+70)/30)^2))+49;
    end
    s_inf=(1/(1+exp((V+(45.3))/(6.8841))));
    s=s+dt*((s_inf-s)/tau_s);
    s_slow_inf=(1/(1+exp((V+(45.3))/(6.8841))));
    s_slow=s_slow+dt*((s_slow_inf-s_slow)/tau_s_slow);
    I_t=(((G_t*r)*((a*s)+(b*s_slow)))*(V-E_K));
    
    % steady_state_outward_K_currents
    r_ss_inf=(1/(1+exp((V+(11.5))/((-1)*(11.82)))));
    tau_r_ss=((10000)/((45.16*exp(0.03577*(V+(50))))+(98.9*exp(((-1)*0.1)*(V+(38))))));
    r_ss=r_ss+dt*((r_ss_inf-r_ss)/tau_r_ss);
    
    s_ss_inf=(1/(1+exp((V+(87.5))/(10.3))));
    tau_s_ss=2100;
    s_ss=s_ss+dt*((s_ss_inf-s_ss)/tau_s_ss);
    
    I_ss=(((G_ss*r_ss)*s_ss)*(V-E_K));
    
    % hyperpolarisation_activated_I_f_current
    y_inf=(1/(1+exp((V+(138.6))/(10.48))));
    tau_y=((1000)/((0.11885*exp((V+(80))/(28.37)))+(0.5623*exp((V+(80))/((-1)*(14.19))))));
    y=y+dt*((y_inf-y)/tau_y);
    
    f_K=(1-f_Na);
    I_f_Na=(((G_f*y)*f_Na)*(V-E_Na));
    I_f_K=(((G_f*y)*f_K)*(V-E_K));
    I_f=(I_f_Na+I_f_K);
    
    % Inward_rectifer_I_k1_current
    I_K1=((((((48)/(exp((V+(37))/(25))+exp((V+(37))/((-1)*(25)))))+...
        (10))*0.001)/(1+exp((V-(E_K+(76.77)))/((-1)*(17)))))+((G_K1...
        *(V-(E_K+(1.73))))/((1+exp(((1.613*F)*(V-(E_K+(1.73))))/...
        (R*T)))*(1+exp((K_o-0.9988)/((-1)*0.124))))));
    
    % sodium_potassium_pump
    sigma=((exp(Na_o/(67.3))-1)/7);
    I_NaK=((((((I_NaK_max*1)/((1+(0.1245*exp(((((-1)*0.1)*V)*F)/(R*T))))+((0.0365*sigma)*exp((((-1)*V)*F)/(R*T)))))*K_o)/(K_o+K_m_K))*1)/(1+((K_m_Na/Na_i)^1.5)));
    
    % Na_Ca_ion_exchanger_current
    I_NaCa=((K_NaCa*((((Na_i^3)*Ca_o)*exp((0.03743*V)*gamma_NaCa))-(((Na_o^3)*Ca_i)*exp((0.03743*V)*(gamma_NaCa-1)))))/(1+(d_NaCa*((Ca_i*(Na_o^3))+(Ca_o*(Na_i^3))))));
    
    % sarcolemmal_calcium_pump_current
    I_pCa=((I_pCa_max*Ca_i)/(Ca_i+(4E-4)));
    
    % L_type_Ca_channelcurrent
    
    d_inf=(1/(1+exp((V+(15.3))/((-1)*(5)))));
    tau_d=((((3.05)*exp(((-1)*0.0045)*((V+(7))^2)))+((1.05)*exp(((-1)*0.002)*((V-(18))^2))))+(0.25));
    d=d+dt*((d_inf-d)/tau_d);
    f_11_inf=(1/(1+exp((V+(26.7))/(5.4))));
    tau_f_11=(((((105)*exp((-1)*(((V+(45))/(12))^2)))+((40)/(1+exp((((-1)*V)+(25))/(25)))))+((15)/(1+exp((V+(75))/(25)))))+1.7);
    f_11=f_11+dt*((f_11_inf-f_11)/tau_f_11);
    f_12_inf=(1/(1+exp((V+(26.7))/(5.4))));
    tau_f_12=(((((41)*exp((-1)*(((V+(47))/(12))^2)))+((80)/(1+exp((V+(55))/((-1)*(5))))))+((15)/(1+exp((V+(75))/(25)))))+(1.7));
    f_12=f_12+dt*((f_12_inf-f_12)/tau_f_12);
    Ca_inact_inf=(1/(1+(Ca_ss/(0.01))));
    Ca_inact=Ca_inact+dt*((Ca_inact_inf-Ca_inact)/tau_Ca_inact);
    I_CaL=(((G_Ca_L*d)*(((0.9+(Ca_inact/10))*f_11)+((0.1-(Ca_inact/10))*f_12)))*(V-65));
    
    %    update concentrationsk
    Na_i = Na_i-dt*(I_Na+I_bNa+I_NaCa*3+I_NaK*3+I_f_Na)*inverseVcF;
    K_i = K_i-dt*(I_stim+I_ss+I_bK+I_t+I_K1+I_f_K-2*I_NaK)*inverseVcF;
    kCa_SR = maxsr-((maxsr-minsr)/(1+(EC/Ca_SR)*(EC/Ca_SR)));
    k1 = k1_/kCa_SR;k2=k2_*kCa_SR;
    dRR = k4*(1-sRR)-k2*Ca_ss*sRR;
    sRR = sRR+dt*dRR;
    sOO = k1*Ca_ss*Ca_ss*sRR/(k3+k1*Ca_ss*Ca_ss);
    Irel = Vrel*sOO*(Ca_SR-Ca_ss);
    Ileak = Vleak*(Ca_SR-Ca_i);
    Iup = Vmaxup/(1.+((Kup*Kup)/(Ca_i*Ca_i)));
    Ixfer = Vxfer*(Ca_ss-Ca_i);
    CaCSQN = Bufsr*Ca_SR/(Ca_SR+Kbufsr);
    dCa_SR = dt*(Iup-Irel-Ileak);
    bjsr = Bufsr-CaCSQN-dCa_SR-Ca_SR+Kbufsr;
    cjsr = Kbufsr*(CaCSQN+dCa_SR+Ca_SR);
    Ca_SR = (sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
    Ca_ssBuf = Bufss*Ca_ss/(Ca_ss+Kbufss);
    dCa_ss = dt*(-Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-I_CaL*inversevssF2));
    bcss = Bufss-Ca_ssBuf-dCa_ss-Ca_ss+Kbufss;
    ccss = Kbufss*(Ca_ssBuf+dCa_ss+Ca_ss);
    Ca_ss = (sqrt(bcss*bcss+4*ccss)-bcss)/2;
    CaBuf = Bufc*Ca_i/(Ca_i+Kbufc);
    dCa_i = dt*((-(I_bCa+I_pCa-2*I_NaCa)*inverseVcF2)-(Iup-Ileak)*(Vsr/Vc)+Ixfer);
    bc = Bufc-CaBuf-dCa_i-Ca_i+Kbufc;
    cc = Kbufc*(CaBuf+dCa_i+Ca_i);
    Ca_i = (sqrt(bc*bc+4*cc)-bc)/2;
    
    %--------------------------------------------------------------------------
    %           Computer active fiblast current
    %--------------------------------------------------------------------------
        fE_Na=Rtonf*log(fNa_o/fNa_i);
    
    % Time- and Voltage- Dependent K+ Current,I_Kv(pA/pF)
    fI_Kv=fG_Kv*fr_Kv*fs_Kv*(v_f-fE_K);
    ft_r=20.3+138*exp(-((v_f+20)/25.9)^2);
    ft_s=1574+5268*exp(-((v_f+23.0)/22.7)^2);
    fr_ioo=1.0/(1.0+exp(-(v_f+20.0)/11));
    fs_ioo=1.0/(1.0+exp((v_f+23.0)/7));
    fr_Kv=fr_Kv+(fr_ioo-fr_Kv)*dt/ft_r;
    fs_Kv=fs_Kv+(fs_ioo-fs_Kv)*dt/ft_s;
    % Timei-Independent Inward-Rectifer K+ Current, I_K1(pA/pF)
    falpha_K1=0.1/(1.0+exp(0.06*(v_f-fE_K-200)));
    fbate_K1=(3.0*exp(0.0002*(v_f-fE_K+100))+exp(0.1*(v_f-fE_K-10)))/...
        (1.0+exp(-0.5*(v_f-fE_K)));
    ref_fIK1=falpha_K1/(falpha_K1+fbate_K1);
    fI_K1=fG_K1*ref_fIK1*(v_f-fE_K);
    
    % Na+ - K+ Exchanger Current, I_NaK(pA/pF)
    fI_NaK=fI_NaK_inf*(fK_o/(fK_o+fK_mk))*((v_f-fV_rev)/(v_f-fB))*...
        ((fNa_i/(fNa_i+fK_mNa))^1.5);
    
    % Background Na+ current, Ib,Na(pA/pF)
    fI_bNa=fG_bNa*(v_f-fE_Na);%
    %--------------------------------------------------------------------------
    %           Computer myocyte-active firboroblast coupling
    %--------------------------------------------------------------------------
    I_tot=(I_Na+I_CaL+I_t+I_ss+I_f+I_K1+I_b+I_NaK+I_NaCa+I_pCa+I_stim)/Cm;
    fI=(fI_Kv+fI_K1+fI_NaK+fI_bNa);
    
    
    
    I_inter_fib=(v_f-V)*G_gap;
    I_inter_myo=Num_Fb*(V-v_f)*G_gap/132;
    if Num_Fb==0
        v_f=v_f-dt*(fI);
    elseif Num_Fb>0
        v_f=v_f-dt*(fI+I_inter_fib/6.3);
    end
    svol=V-dt* (I_tot+I_inter_myo);
    
    if(time > (stims1*(xns-2)))&& (time < (stims1*(xns-1)+tbegin))
        APV = [APV; V];
    end
    %--------------------------------------------------------------------------
    %           Computer APD
    %--------------------------------------------------------------------------
    if(time > (stims1*(xns-2)))&& (time < (stims1*(xns-1)+tbegin))
        APV = [APV; V];
    end
    
    if(time > (stims1*(xns-1)+tbegin))&& (time < (stims1*(xns-1)+sdur+tbegin))
        AP90=max(APV)-(max(APV)-min(APV))*0.9;
    end
    
    if(time > (stims1*(xns-1)+sdur+tbegin))
        if ((svol-AP90)*(V-AP90)<0)
            napdd=napdd+1;
            apdtime(napdd) = time;
        end
        if(napdd==3)
            DI=apdtime(2)-apdtime(1);
            APD=apdtime(3)-apdtime(2);
        end
    end
    V = svol;
    
    %     if time>=25
    if mod(step,100)==0
        st=[st,time];
        sv=[sv, V];
        svf =[svf, v_f];
        sCa=[sCa, Ca_i*1000];
        sa=[sa, Na_i];
        sb=[sb,K_i];
    end
    %     end
    step=step+1;
    time=time+dt;
end
if(napdd==3)
    PD = [PD;APD];
    PI =[PI; DI];
    cycl =[cycl; xns];
end
disp(['APD90: ',num2str(PD)]);
figure
subplot(2,1,1)
hold on
plot(st,sv)
title(' Cardiomyocyte V_m (mV)','Fontsize',15);
subplot(2,1,2)
hold on
plot(st,svf)
title(' Active fibroblast V_f (mV)','Fontsize',15);
figure
subplot (3,1,1)
hold on
plot(st,sCa)
title(' Ca_i(nM)','Fontsize',15);
subplot (3,1,2)
plot(st,sa)
title('Na_i(uM)','Fontsize',25);
subplot (3,1,3)
plot(st,sb)
title('K_i(uM)','Fontsize',15);
AA1 =[st; sv;svf;sCa]';
dlmwrite('Pandit_MacCannell.txt',AA1,'delimiter','\t','precision',6)
t2 = clock;
disp(['total time: ',num2str(etime(t2,t1))]);
% dlmwrite('ssPan20150908.txt',M1,'delimiter','\t','precision',6)