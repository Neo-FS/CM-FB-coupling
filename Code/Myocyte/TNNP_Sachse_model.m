% MODEL NAME: TNNP_Sachse in humans
% SHORT DESCRIPTION: This model reproduces the action potential recorded experimentally
% for myocytes isolated from the adult human left ventricle (Ten Tusscher et al. 2004) and
% active fibroblasts from Sachse (Sachse et al. 2008).
% Maximum INaK current,I_NaK_max, and Scaling factor for INaCa , K_NaCa, are adjusted to
% satisfy the convergencein of Na+，K+ and Ca2+ in epicardial myocytes, M myocytes and endocardial myocytes.

clc;
t1=clock;
%------------------------------------------------------------------------
%                  Parameter for human ventricular myocytes
%------------------------------------------------------------------------
Num_Fb = 0;  % the numeb of active fibroblast
Celltype = 1; % 1,Epicardial cell; 2,Endocardial cell;  3, Middle
%  Physical constant
R=8314.5;  % millijoule_per_mole_kelvin (in membrane) - Ideal gas constant
T=310;     % kelvin (in membrane) - Absolute tempehumanure
F=96487;   % coulomb_per_mole (in membrane) - Faraday constant
Rtonf = (R*T)/F;

% Cellular capacitance
C_m=0.185; % microF (in membrane) - human myocyte membrane capacitance

% External concenthumanions
K_o=5.4; % millimoles per liter
Ca_o=2.0; % millimoles per liter
Na_o=140.0; % millimoles per liter

% Initial internal concenthumanions
Ca_i=0.00008; % millimoles per liter
Ca_sr=0.56; % millimoles per liter
Na_i=11.6; % millimoles per liter
K_i=138.3; % millimoles per liter

%Intracellular volumes
V_c=0.016404; %  nanoliliter - Cellular volumes of human ventricular cells
V_sr=0.001094; %  nanoliliter - endoplasmic reticulum volumes of human ventricular cells

%Parameters for currents
G_Na=14.838; % nanoS per picoF;
G_bNa=0.00029; % nanoS per picoF;
G_CaL=0.000175; % nanoS per picoF;
G_bCa=0.000592; % nanoS per picoF;
G_Kr=0.096; % nanoS per picoF;
p_KNa=0.03;
G_K1=5.405; % nanoS per picoF;
if Celltype == 1
    G_Ks=0.245; % nanoS per picoF;
    G_to=0.294; % nanoS per picoF;
    P_NaK=0.954*1.362;  % picoS per  picoF;
    K_NaCa=0.945*1000.0; % picoS per  picoF;
elseif Celltype == 2
    G_Ks=0.245; % nanoS per picoF;
    G_to=0.073; % nanoS per picoF;
    P_NaK=0.886*1.362;  % picoS per  picoF;
    K_NaCa=0.880*1000.0; % picoS per  picoF;
elseif Celltype == 3
    G_Ks=0.062; % nanoS per picoF;
    G_to=0.294; % nanoS per picoF;
    P_NaK=0.966*1.362;  % picoS per  picoF;
    K_NaCa=1.320*1000.0; % picoS per  picoF;
else
    error('Error:Incorrect input of parameters')
end
K_mK=1.0; % millimoles per liter
K_mNa=40.0; % millimoles per liter
K_mNai=87.5; % millimoles per liter
K_mCa=1.38; % millimoles per liter
k_sat=0.1;
gama=0.35;
G_pCa=0.825; % nanoS per picoF;
K_pCa=0.0005; % nanoS per picoF;
G_pK=0.0146; % nanoS per picoF;
alpha=2.5;

%Calcium dynamics
Buf_c=0.15;  % millimoles per liter
K_bufc=0.001;  % millimoles per liter
Buf_sr=10.0;  % millimoles per liter
K_bufsr=0.3;  % millimoles per liter
V_maxup=0.000425;  % millimoles per liter per millisecond
K_up=0.00025; % millimoles per liter

% Initial internal variable of myocytes
v=-82.125080404886690;
m=0.003299066621003;
h=0.649393306465186;
j=0.648624810520527;
x_r1=3.429794493560507e-04;
x_r2=0.439106523501946;
x_s=0.004034096457745;
r=0;
s=1;
d=3.419684265587707e-05;
f=0.999774079145077;
f_Ca=0.321450251716759;
g=0.999970129256348;

%--------------------------------------------------------------------------
%              Parameter for Sachse active fibroblast model
%--------------------------------------------------------------------------
fC_m=4.5; % picroF (in membrane) - active firboblast  membrane capacitance
%  Physical constants
fF= 96500; fR=8310;  fT=295;
fK_o=5;  % millimoles per liter
fNa_o=Na_o; % millimoles per liter
fK_i= 140; % millimoles per liter
fNa_i=10; % millimoles per liter

G_gap=3.0;  % nanoS
G_Kir=1.02; %  nanoS per picoF;
a_Kir=0.94;
b_Kir=1.26;

P_Shkr=5.4e-6;
kvo=30e-3;%
k_vo=2e-3;%
zv=1.28;%
z_v=-1.53;%
ko=77e-3;%
k_o=18e-3;%
fG_b=6.9e-3;
fE_b=0.0;


fE_K=-87;  % milliVolt
fV_rev=-150.0; % milliVolt
fB=-200; % milliVolt
fK_mk=1.0;% millimoles per liter
fK_mNa=11.0;% millimoles per liter
fI_NaK_inf=2.002; % picoA per picoF;

% Initial internal variable of active fibroblast
v_f=-49.6;
fr_Kv=0;
fs_Kv=1;

%  Initial Values of fibroblast model for [K_o]=5mM
V_f = -59.5021795;%mV
C_0Shkr = 0.920473725503476;
C_1Shkr = 0.076853417684277;
C_2Shkr = 0.002406280448088;
C_3Shkr = 3.348472843107640e-05;
C_4Shkr = 1.747344724783623e-07;
O_Shkr = 7.474752795471331e-07;

%--------------------------------------------------------------------------
%             Parameter for simulation duration
%--------------------------------------------------------------------------
PD =[];  PI =[];cycle =[];
% for xns = 1:15
xns =10;
dt = 0.005;
sdur = 1.0;
I_stim = 0;
stimstrength =-52;  %  stimstrength
stims1 = 1000;
stims2 = 1000;
tbegin=20;
xnstims1 = stims1+tbegin+sdur;%+tbegin;
xnstims2 =stims1*(xns-1)+tbegin+stims2+sdur;%tbegin;
endtime=xnstims2+stims2;  %duhumanion of the simulation
nswitch = 0;
ncounts1=1;
napdd=0;
apdtime = zeros(xns,1);
time=0; step=0; %timestep(ms)
st=[];sv=[];sCa = []; sfv = [];sa=[];sb=[];sc=[];sd=[];se=[];sf=[];sg=[];APV = [];
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
    %           Computer myocytes current
    %--------------------------------------------------------------------------
    %  reversal potentials
    E_K=Rtonf*(log(K_o/K_i));
    E_Na=Rtonf*(log(Na_o/Na_i));
    E_Ks=Rtonf*(log((K_o+p_KNa*Na_o)/(K_i+p_KNa*Na_i)));
    E_Ca=0.5*Rtonf*(log(Ca_o/Ca_i));
    alpha_K1=0.1/(1.0+exp(0.06*(v-E_K-200)));
    bata_K1=(3.0*exp(0.0002*(v-E_K+100))+...
        exp(0.1*(v-E_K-10)))/(1.0+...
        exp(-0.5*(v-E_K)));
    
    rec_I_NaK=(1.0/(1.0+0.1245*exp(-0.1*v/Rtonf)+0.0353*exp(-v/Rtonf)));
    rec_I_pK=1.0/(1.0+exp((25-v)/5.98));
    
    %  I_Na_current
    alpha_m=1.0/(1.0+exp((-60.0-v)/5.0));
    bate_m=0.1/(1+exp((v+35.0)/5.0))+0.1/(1.0+exp((v-50.0)/200.0));
    t_m=alpha_m*bate_m;
    m_inf=1.0/((1.0+exp((-56.86-v)/9.03))^2);
    if v>=-40.0
        alpha_h=0.0;
        bate_h=0.77/(0.13*(1.0+exp(-(v+10.66)/11.1)));
        alpha_j=0.0;
        bate_j=0.6*exp(0.057*v)/(1.0+exp(-0.1*(v+32.0)));
    else
        alpha_h=0.057*exp(-(v+80.0)/6.8);
        bate_h=2.7*exp(0.079*v)+3.1e5*exp(0.3485*v);
        alpha_j=(-2.5428e4*exp(0.2444*v)-(6.948e-6)*exp(-0.04391*v))*...
            (v+37.78)/(1+exp(0.311*(v+79.23)));
        bate_j=0.02424*exp(-0.01052*v)/(1.0+exp(-0.1378*(v+40.14)));
    end
    t_h=1.0/(alpha_h+bate_h);
    t_j=1.0/(alpha_j+bate_j);
    h_inf=1.0/((1.0+exp((v+71.55)/7.43))^2);
    j_inf=h_inf;
    m=m_inf-(m_inf-m)*exp(-dt/t_m);
    h=h_inf-(h_inf-h)*exp(-dt/t_h);
    j=j_inf-(j_inf-j)*exp(-dt/t_j);
    I_Na=G_Na*(m^3)*h*j*(v-E_Na);
    
    %  I_CaL_current
    d_inf=1.0/(1.0+exp((-5.0-v)/7.5));
    alpha_d=1.4/(1.0+exp((-35.0-v)/13.0))+0.25;
    bate_d=1.4/(1+exp((5.0+v)/5.0));
    gama_d=1.0/(1.0+exp((50.0-v)/20.0));
    t_d=alpha_d*bate_d+gama_d;
    
    f_inf=1.0/(1.0+exp((v+20.0)/7.0));
    t_f=1125.0*exp(-((v+27.0)^2)/300)+80.0+165.0/(1+exp((25.0-v)/10.0));
    
    alpha_fCa=1.0/(1.0+(Ca_i/0.000325)^8);
    bate_fCa=0.1/(1.0+exp((Ca_i-0.0005)/0.0001));
    gama_fCa=0.2/(1.0+exp((Ca_i-0.00075)/0.0008));
    f_Ca_inf=(alpha_fCa.*bate_fCa+gama_fCa+0.23)/1.46;
    if Ca_i<=0.00035
        g_inf=1/(1+(Ca_i/0.00035)^6);
    else
        g_inf=1/(1+(Ca_i/0.00035)^16);
    end
    d=d_inf-(d_inf-d)*exp(-dt/t_d);
    f=f_inf-(f_inf-f)*exp(-dt/t_f);
    f_Caold=f_Ca;
    f_Ca=f_Ca_inf-(f_Ca_inf-f_Ca)*exp(-dt/2);
    if (f_Ca>f_Caold)&&(v>-37)
        f_Ca=f_Caold;
    end
    g_old=g;
    g=g_inf-(g_inf-g)*exp(-dt/2);
    if g>g_old&&v>-37
        g=g_old;
    end
    I_CaL=G_CaL*d*f*f_Ca*4*v*(F/Rtonf)*(Ca_i*exp(2*v/Rtonf)-0.341*Ca_o)...
        /(exp(2*v/Rtonf)-1.0);
    
    %  I_to_current
    r_inf=1.0/(1.0+exp((20.0-v)/6.0));
    t_r=9.5*exp((-(v+40.0)^2)/1800.0)+0.8;
    if Celltype == 2
        s_inf=1.0/(1.0+exp((20.0+v)/5.0));
        t_s=85.0*exp((-(v+45.0)^2)/320.0)+5.0/(1.0+exp((v-20.0)/5.0))+3.0;
    elseif Celltype == 1|Celltype == 3
        s_inf=1/(1+exp((28+v)/5));
        t_s=1000*exp(-(v+67)^2/1000)+8;
    end
    
    s=s_inf-(s_inf-s)*exp(-dt/t_s);
    r=r_inf-(r_inf-r)*exp(-dt/t_r);
    I_to=G_to*r*s*(v-E_K);
    
    %  I_Kr_current
    x_r1_inf=1.0/(1.0+exp((-26.0-v)/7.0));
    alpha_xr1=450.0/(1.0+exp((-45.0-v)/10.0));
    bate_xr1=6.0/(1.0+exp((v+30.0)/11.5));
    t_xr1=alpha_xr1*bate_xr1;
    x_r2_inf=1.0/(1.0+exp((v+88.0)/24.0));
    alpha_xr2=3.0/(1+exp((-60.0-v)/20.0));
    bate_xr2=1.12/(1+exp((v-60.0)/20.0));
    t_xr2=alpha_xr2*bate_xr2;
    x_r1=x_r1_inf-(x_r1_inf-x_r1)*exp(-dt/t_xr1);
    x_r2=x_r2_inf-(x_r2_inf-x_r2)*exp(-dt/t_xr2);
    I_Kr = G_Kr*(sqrt(K_o/5.4))*x_r1*x_r2*(v-E_K);
    
    %  I_Ks_current
    x_s_inf=1.0/(1+exp((-5.0-v)/14.0));
    alpha_xs=1100.0/(sqrt(1.0+exp((-10.0-v)/6.0)));
    bate_xs=1.0/(1.0+exp((v-60.0)/20.0));
    t_xs = alpha_xs*bate_xs;
    x_s=x_s_inf-(x_s_inf-x_s)*exp(-dt/t_xs);
    I_Ks = G_Ks*((x_s)^2)*(v-E_Ks);
    
    
    x_K1_inf = alpha_K1/(alpha_K1+bata_K1);
    I_K1 = G_K1*x_K1_inf*(v-E_K);
    
    I_NaCa = K_NaCa*(1.0/((K_mNai)^3+(Na_o)^3))*(1.0/(K_mCa+Ca_o))*...
        (1.0/(1+k_sat*exp((gama-1)*v/Rtonf)))*...
        (exp(gama*v/Rtonf)*(Na_i)^3*Ca_o-...
        exp((gama-1)*v/Rtonf)*(Na_o)^3*Ca_i*alpha);
    
    %  I_K1_current
    I_NaK= P_NaK*(K_o/(K_o+K_mK))*(Na_i/(Na_i+K_mNa))*rec_I_NaK;
    
    %  I_pCa_current
    I_pCa=G_pCa*Ca_i/(K_pCa+Ca_i);
    
    %  I_pK_current
    I_pK=G_pK*rec_I_pK*(v-E_K);
    
    %  I_b_current
    I_bNa=G_bNa*(v-E_Na);
    I_bCa=G_bCa*(v-E_Ca);
    
    % Intracellular dynamic ion concentration
    I_Ca=-((I_CaL+I_bCa+I_pCa-2*I_NaCa)*C_m)/(2*V_c*F);
    I_rel=((0.016464*Ca_sr*Ca_sr)/(0.0625+Ca_sr*Ca_sr)+0.008232)*d*g;
    I_leak=0.00008*(Ca_sr-Ca_i);
    I_up=V_maxup/(1.0+(K_up*K_up)/(Ca_i*Ca_i));
    I_Casr=I_up-I_rel-I_leak;
    
    Ca_srbufsr=Buf_sr*Ca_sr/(Ca_sr+K_bufsr);
    dCa_sr=dt*(V_c/V_sr)*I_Casr;
    bjsr=Buf_sr-Ca_srbufsr-dCa_sr-Ca_sr+K_bufsr;
    cjsr=K_bufsr*(Ca_srbufsr+dCa_sr+Ca_sr);
    Ca_sr=(sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
    Ca_ibufc=Buf_c*Ca_i/(Ca_i+K_bufc);
    dCa_i=dt*(I_Ca-I_Casr);
    bc=Buf_c-Ca_ibufc-dCa_i-Ca_i+K_bufc;
    cc=K_bufc*(Ca_ibufc+dCa_i+Ca_i);
    Ca_i=(sqrt(bc*bc+4*cc)-bc)/2;
    
    Na_i=Na_i-dt*(I_Na+I_bNa+3*I_NaK+3*I_NaCa)*C_m/(V_c*F);
    K_i=K_i-dt*(I_stim+I_K1+I_to+I_Kr+I_Ks-2*I_NaK+I_pK)*C_m/(V_c*F);
    
    %--------------------------------------------------------------------------
    %           Computer active fiblast current
    %--------------------------------------------------------------------------
    %  I_Kir current
    fE_K = Rtonf*log(fK_o/fK_i);
    O_Kir=1/(a_Kir+exp(b_Kir*(V_f-fE_K)/Rtonf));
    fI_Kir=G_Kir*O_Kir*sqrt(fK_o*0.001)*(V_f-fE_K);
    
    %  I_Shkr current
    kv=kvo*exp(V_f*zv/Rtonf);
    k_v=k_vo*exp(V_f*z_v/Rtonf);
    if C_0Shkr<0
        C_0Shkr=0;
    elseif C_0Shkr>1
        C_0Shkr=1;
    else
        dC_0Shkr=(dt)*(k_v*C_1Shkr-4*kv*C_0Shkr);
        C_0Shkr=C_0Shkr+dC_0Shkr;
    end
    if C_1Shkr<0
        C_1Shkr=0;
    elseif C_1Shkr>1
        C_1Shkr=1;
    else
        dC_1Shkr=(dt)*(2*k_v*C_2Shkr+4*kv*C_0Shkr-(3*kv+k_v)*C_1Shkr);
        C_1Shkr=C_1Shkr+dC_1Shkr;
    end
    if C_2Shkr<0
        C_2Shkr=0;
    elseif C_2Shkr>1
        C_2Shkr=1;
    else
        dC_2Shkr=(dt)*(3*k_v*C_3Shkr+3*kv*C_1Shkr-(2*kv+2*k_v)*C_2Shkr);
        C_2Shkr=C_2Shkr+dC_2Shkr;
    end
    if C_3Shkr<=0
        C_3Shkr=0;
    elseif C_3Shkr>=1
        C_3Shkr=1;
    else
        dC_3Shkr=(dt)*(4*k_v*C_4Shkr+2*kv*C_2Shkr-(kv+3*k_v)*C_3Shkr);
        C_3Shkr=C_3Shkr+dC_3Shkr;
    end
    if C_4Shkr<=0
        C_4Shkr=0;
    elseif C_4Shkr>=1
        C_4Shkr=1;
    else
        dC_4Shkr=(dt)*(k_o*O_Shkr+kv*C_3Shkr-(ko+4*k_v)*C_4Shkr);
        C_4Shkr=C_4Shkr+dC_4Shkr;
    end
    if O_Shkr<0
        O_Shkr=0;
    elseif O_Shkr>1
        O_Shkr=1;
    else
        dO_Shkr=(dt)*(ko*C_4Shkr-k_o*O_Shkr);
        O_Shkr=O_Shkr+dO_Shkr;
    end
    fI_Shkr=P_Shkr*O_Shkr*(V_f*fF/Rtonf)*(fK_i-fK_o*exp(-V_f/Rtonf))/...
        (1-exp((-1*V_f)/Rtonf));
    
    %  I_b current
    fI_b=fG_b*(V_f-fE_b);
    
    %--------------------------------------------------------------------------
    %           Computer myocyte-active firboroblast coupling
    %--------------------------------------------------------------------------
    %Total Membrance Current and Memberance Potential(pA/pF):
    I_tot=I_Kr+I_Ks+I_K1+I_to+I_Na+I_bNa+I_CaL+I_bCa+I_NaCa+I_NaK+I_pCa+...
        I_pK+I_stim;
    fI_tot=fI_Kir+fI_Shkr+fI_b;
    
    % coupling current
    I_gap=G_gap*(V_f-v);
    
    % coupling membrane potential
    if Num_Fb ==0
        V_f=V_f-(dt)*(fI_tot)/fC_m;
    elseif Num_Fb>0
        V_f=V_f-(dt)*(fI_tot+I_inter_fib)/fC_m;
    end
    V_f=V_f-dt*(fI_tot+I_gap/fC_m);
    svol=v-dt*(I_tot-Num_Fb*I_gap/185);%%%%C_m=72;
    
    % computer APD
    if(time > (stims1*(xns-2)))&& (time < (stims1*(xns-1)+tbegin))
        APV = [APV; svol];
    end
    if(time > (stims1*(xns-1)+tbegin))&& (time < (stims1*(xns-1)+sdur+tbegin))
        AP90=max(APV)-(max(APV)-min(APV))*0.9;
    end
    if(time > (stims1*(xns-1)+sdur+tbegin))
        if ((svol-AP90)*(v-AP90)<0)
            napdd=napdd+1;
            apdtime(napdd) = time;
        end
        if(napdd==3)
            DI=apdtime(2)-apdtime(1);
            APD=apdtime(3)-apdtime(2);
        end
    end
    v = svol;
    
    %--------------------------------------------------------------------------
    %                              Write the PARAMEITER
    %--------------------------------------------------------------------------
    if mod(step,100)==0
        st=[st,time];
        sv=[sv,v];
        sfv = [sfv, v_f];
        sCa = [sCa, Ca_i*1000];
        sa=[sa, Na_i];
        sb=[sb, K_i];
    end
    step=step+1;
    time = time +dt;
end
if(napdd==3)
    PD =[PD;APD];
    PI =[PI; DI];
    cycle =[cycle; xns];
end
disp(['APD90: ',num2str(PD)]);
figure
subplot(2,1,1)
plot(st,sv)
title(' Cardiomyocyte V_m (mV)','Fontsize',15);
subplot(2,1,2)
plot(st,sfv)
title(' Active fibroblast V_f (mV)','Fontsize',15);

figure
subplot(3,1,1)
plot(st,sCa)
title(' Ca_i(uM)','Fontsize',15);
subplot(3,1,2)
plot(st,sa)
title('Na_i(uM)','Fontsize',15);
subplot(3,1,3)
plot(st,sb)
title('K_i(uM)','Fontsize',15);
AA1=[st;sv;sfv]';
dlmwrite('TNNP_Sachse.txt',AA1,'delimiter','\t','precision',6)
t2 = clock;
disp(['total time: ',num2str(etime(t2,t1))]);