global RHOV_sur RHOV_A Resis_a Resis_s P_Va Velo_fric Theta_LL_sur  % RHOV_sur and Theta_L_sur should be stored at each time step.
global z_ref z_srT z_srm VK_Const d0_disp U_wind MO Ta U Ts Zeta_MO Stab_m Stab_T       % U_wind is the mean wind speed at height z_ref (m��s^-1), U is the wind speed at each time step.
global Rv g HR_a NL NN Evap KT RHOV Theta_LL EVAP Evapo
global Tao LAI Tp_t Trap
global H1 H2 H3 H4 alpha_h bx LR Lm Srt Elmn_Lnth DeltZ RL TIME rwuef hh
global PT_PM_VEG PE_PM_SOIL
global r_s_VEG r_s_SOIL r_a_VEG r_a_SOIL Rn_SOIL Rn w e0_Ts e_a_Ts e0_Ta e_a wt wc Srt_1 rl LAI_act rl_min
global hh_v T TT Theta_UU_sur Theta_UU alfa_s JN coszen1 Coef_uV ALBD ALBD_V ALBD_V_VIS ALBD_VB_vis ALBD_VD_vis ALBD_G Albd_g_vis Albd_g_nir ALBD_V_NIR ALBD_VB_nir ALBD_VD_nir

%%%%%%% LAI and light extinction coefficient calculation %%%%%%%%%%%%%%%%%%
LAI(KT)=1.2;
if LAI(KT)<=2
    LAI_act(KT)=LAI(KT);
elseif LAI(KT)<=4
    LAI_act(KT)=2;
else
    LAI_act(KT)=0.5*LAI(KT);
end
LAI_act(KT)=0.23;
Tao=0.56;  %light attenual coefficient
% Set constants
sigma = 4.903e-9; % Stefan Boltzmann constant MJ.m-2.day-1 FAO56 pag 74
lambdav = 2.45;    % latent heat of evaporation [MJ.kg-1] FAO56 pag 31
% Gieske 2003 pag 74 Eq33/DKTgman 2002
% lambda=2.501-2.361E-3*t, with t temperature evaporative surface (?C)
% see script Lambda_function_t.py
Gsc = 0.082;      % solar constant [MJ.m-2.mKT-1] FAO56 pag 47 Eq28
eps = 0.622;       % ratio molecular weigth of vapour/dry air FAO56 p26 BOX6
R = 0.287;         % specific gas [kJ.kg-1.K-1]    FAO56 p26 box6
Cp = 1.013E-3;     % specific heat at cte pressure [MJ.kg-1.?C-1] FAO56 p26 box6
k = 0.41;          % karman's cte   []  FAO 56 Eq4
Z=3421;             % altitute of the location(m)
as=0.25;           % regression constant, expressKTg the fraction of extraterrestrial radiation FAO56 pag50
bs=0.5;
alfa=0.23;         % albeo of vegetation set as 0.23
z_m=10;            % observation height of wKTd speed; 10m
Lz=240*pi()/180;   % latitude of Beijing time zone west of Greenwich
Lm=252*pi()/180;    % latitude of Local time, west of Greenwich
% albedo of soil calculation;
Theta_LL_sur(KT)=Theta_LL(NL,1);
Theta_UU_sur(KT)=Theta_UU(NL,1);
if Theta_LL_sur(KT)<0.1
    alfa_s(KT)=0.25;
elseif Theta_LL_sur(KT)<0.25
    alfa_s(KT)=0.35-Theta_LL_sur(KT);
else
    alfa_s(KT)=0.1;
end
%% Albedo of land surface

YEAR=2015;
MONTH=12;
DAY=1;
HOUR=9.5;

JN(KT)=fix((TIME/3600+HOUR)/24)+335;    % day number
if JN(KT)>365
    JN(KT)=JN(KT)-365;
end

%% RADIATION PARAMETERS CALCULATION
% compute dr - KTverse distance to the sun
% [rad]
% FAO56 pag47 Eq23
dr(KT) = 1+0.033*cos(2*pi()*JN(KT)/365);

% compute delta - solar declKTation
% [rad]
% FAO56 pag47 Eq24
delta(KT) = 0.409*sin(2*pi()*JN(KT)/365-1.39);

% compute Sc - seasonnal correction of solar time
% [hour]
% FAO56 pag47 Eq32
Sc = [];
b(KT) = 2.0*pi()*(JN(KT)-81.0)/364.0;    % Eq 34
Sc(KT) = 0.1645*sin(2*b(KT)) - 0.1255*cos(b(KT)) - 0.025*sin(b(KT));

% compute w - solar time angle at the midpoKTt of the period (time)
% [rad]
% FAO56 pag48 Eq31
w(KT)=pi()/12*((TIME/3600-fix(TIME/3600/24-0.001)*24-0.5+0.06667*(Lz-Lm)+Sc(KT))-12);
[HRI,COSZEN]=cosen(JN(KT),HOUR,TIME/3600,0,0,34);
if KT==1
    Coef_u=COSZEN;
else
    Coef_u=coszen1(KT-1);%sin(0.599)*sin(delta(KT)) + cos(0.599)*cos(delta(KT))*cos(w(KT));
end
Coef_uV(KT)=Coef_u;%sin(0.599)*sin(delta(KT)) + cos(0.599)*cos(delta(KT))*cos(w(KT));
%%
ALBD_S_vis=[0.15 0.11 0.10 0.09 0.08 0.07 0.06 0.05 0.083];
ALBD_S_nir=[0.30 0.22 0.20 0.18 0.16 0.14 0.12 0.10 0.166];
ALBD_d_vis=[0.27 0.22 0.20 0.18 0.16 0.14 0.12 0.10 0.166];
ALBD_d_nir=ALBD_d_vis.*2;%[0.30 0.22 0.20 0.18 0.16 0.14 0.12 0.10];
Delt_albd=max(0.01*(11-40*Theta_LL_sur(KT)),0);
TYPE_S=9;
Albd_g_vis(KT)=ALBD_S_vis(TYPE_S)+min(Delt_albd,ALBD_S_vis(TYPE_S));
Albd_g_nir(KT)=Albd_g_vis(KT)*2;
ALBD_G(KT)=(Albd_g_vis(KT)+Albd_g_nir(KT))/2;
%% albedo for vegetation clmdoc P47
%for grassland, 10
ALBD_VBF_vis=[0.04	0.04	0.05	0.07	0.06	0.07	0.14	0.06	0.07	0.07	0.06	0.06	0.06	0.06	0.95	0.19];
ALBD_VBF_nir=[0.2	0.2	0.23	0.24	0.24	0.26	0.32	0.21	0.26	0.25	0.18	0.24	0.22	0.22	0.65	0.38];
TYPE_veg=10;
Coef_w_vis=0.15;Coef_w_nir=0.85;
beta=0.5;
SAI=LAI(KT);
ALBD_GB_VIS=Albd_g_vis(KT);
ALBD_GB_NIR=Albd_g_nir(KT);
%     Fc(KT)=1-exp(-1*(Tao*LAI(KT)));
ALBD_VB_vis(KT)=ALBD_VBF_vis(TYPE_veg)*(1-exp(-1*(Coef_w_vis*beta*SAI)/(Coef_u*ALBD_VBF_vis(TYPE_veg))))+ALBD_GB_VIS*exp(-1*(1+0.5/Coef_u)*SAI);
ALBD_VD_vis(KT)=ALBD_VBF_vis(TYPE_veg)*(1-exp(-1*(Coef_w_vis*beta*SAI)/(Coef_u*ALBD_VBF_vis(TYPE_veg))))+ALBD_GB_VIS*exp(-1*(1+0.5/Coef_u)*SAI);
ALBD_VB_nir(KT)=ALBD_VBF_nir(TYPE_veg)*(1-exp(-1*(Coef_w_nir*beta*SAI)/(0.5*ALBD_VBF_nir(TYPE_veg))))+ALBD_GB_NIR*exp(-1*(1+0.5/0.5)*SAI);
ALBD_VD_nir(KT)=ALBD_VBF_nir(TYPE_veg)*(1-exp(-1*(Coef_w_nir*beta*SAI)/(0.5*ALBD_VBF_nir(TYPE_veg))))+ALBD_GB_NIR*exp(-1*(1+0.5/0.5)*SAI);
ALBD_V_VIS(KT)=ALBD_VB_vis(KT)+ ALBD_VD_vis(KT);
ALBD_V_NIR(KT)=ALBD_VB_nir(KT)+ ALBD_VD_nir(KT);
ALBD_V(KT)=(ALBD_V_VIS(KT)+ALBD_V_NIR(KT))/2;
ALBD(KT)=ALBD_G(KT)*(1-0.9)+ALBD_V(KT)*(0.9);
%     ALBD=ALBD_G.*(1-0.69)+ALBD_V.*(0.69);
%% AIR PARAMETERS CALCULATION
% compute DELTA - SLOPE OF SATURATION VAPOUR PRESSURE CURVE
% [kPa.?C-1]
% FAO56 pag 37 Eq13
DELTA(KT) = 4098*(0.6108*exp(17.27*Ta(KT)/(Ta(KT)+237.3)))/(Ta(KT)+237.3)^2;
% ro_a - MEAN AIR DENSITY AT CTE PRESSURE
% [kg.m-3]
% FAO56 pag26 box6
Pa=101.3*((293-0.0065*Z)/293)^5.26;
ro_a(KT) = Pa/(R*1.01*(Ta(KT)+273.16));
% compute e0_Ta - saturation vapour pressure at actual air temperature
% [kPa]
% FAO56 pag36 Eq11
e0_Ta(KT) = 0.6108*exp(17.27*Ta(KT)/(Ta(KT)+237.3));
e0_Ts(KT) = 0.6108*exp(17.27*Ts(KT)/(Ts(KT)+237.3));
% compute e_a - ACTUAL VAPOUR PRESSURE
% [kPa]
% FAO56 pag74 Eq54
e_a(KT) = e0_Ta(KT)*HR_a(KT);
e_a_Ts(KT) = e0_Ts(KT)*HR_a(KT);

% gama - PSYCHROMETRIC CONSTANT
% [kPa.?C-1]
% FAO56 pag31 eq8
gama = 0.664742*1e-3*Pa;

hh_v(KT)=0.12;
rl_min(KT)=100;
%     Rn(KT) =Rn_msr(KT);
if TIME<=1800*3600
    Rn_SOIL(KT) =Rn(KT)*0.68;  % exp(-1*(Tao*LAI(KT)))net radiation for soil
else
    Rn_SOIL(KT) =Rn(KT)*0.68;
end
%% SURFACE RESISTANCE PARAMETERS CALCULATION
R_a=0.81;R_b=0.004*24*11.6;R_c=0.05;
rl(KT)=rl_min(KT)/((R_b*Rn(KT)+R_c)/(R_a*(R_b*Rn(KT)+1)));

% r_s - SURFACE RESISTANCE
% [s.m-1]
% VEG: Dingman pag 208 (canopy conductance) (equivalent to FAO56 pag21 Eq5)
r_s_VEG(KT) = rl(KT)/LAI_act(KT);

% SOIL: equation 20 of van de Griend and Owe, 1994

if KT<=1047
    r_s_SOIL(KT)=10.0*exp(0.3563*100.0*(0.18-Theta_LL_sur(KT)));   % set as minmum soil moisture for potential evaporation
elseif KT<=3624
    r_s_SOIL(KT)=10.0*exp(0.3563*100.0*(0.18-Theta_LL_sur(KT)));   % set as minmum soil moisture for potential evaporation
else
    r_s_SOIL(KT)=10.0*exp(0.3563*100.0*(0.18-Theta_LL_sur(KT)));   % set as minmum soil moisture for potential evaporation
end
% r_a - AERODYNAMIC RESISTANCE
% [s.m-1]
% FAO56 pag20 eq4- (d - zero displacement plane, z_0m - roughness length momentum transfer, z_0h - roughness length heat and vapour transfer, [m], FAO56 pag21 BOX4
r_a_VEG(KT) = log((2-(2*hh_v(KT)/3))/(0.123*hh_v(KT)))*log((2-(2*hh_v(KT)/3))/(0.0123*hh_v(KT)))/((k^2)*U(KT))*100;  % s m-1
% r_a of SOIL
% Liu www.hydrol-earth-syst-sci.net/11/769/2007/
% equation for neutral conditions (eq. 9)
% only function of ws, it is assumed that roughness are the same for any type of soil

RHOV_sur(KT)=RHOV(NN);
P_Va(KT)=0.611*exp(17.27*Ta(KT)/(Ta(KT)+237.3))*HR_a(KT);  %The atmospheric vapor pressure (KPa)  (1000Pa=1000.1000.g.100^-1.cm^-1.s^-2)

RHOV_A(KT)=P_Va(KT)*1e4/(Rv*(Ta(KT)+273.15));              %  g.cm^-3;  Rv-cm^2.s^-2.Cels^-1

z_ref=200;          % cm The reference height of tempearature measurement (usually 2 m)
d0_disp=0;          % cm The zero-plane displacement (=0 m)
z_srT=0.1;          % cm The surface roughness for the heat flux (=0.001m) 0.01m
VK_Const=0.41;   % The von Karman constant (=0.41)
z_srm=0.1;          % cm The surface roughness for momentum flux (=0.001m) 0.01m
U_wind=198.4597; % The mean wind speed at reference height.(cm.s^-1)

MO(KT)=(Ta(KT)*U(KT)^2)/(g*(Ta(KT)-T(NN))*log(z_ref/z_srm));  % Wind speed should be in cm.s^-1, MO-cm;

Zeta_MO(KT)=z_ref/MO(KT);

if abs(Ta(KT)-T(NN))<=0.01
    Stab_m(KT)=0;
    Stab_T(KT)=0;
elseif Zeta_MO(KT)<0
    Stab_T(KT)=-2*log((1+sqrt(1-16*Zeta_MO(KT)))/2);
    Stab_m(KT)=-2*log((1+(1-16*Zeta_MO(KT))^0.25)/2)+Stab_T(KT)/2+2*atan((1-16*Zeta_MO(KT))^0.25)-pi/2;
else
    if Zeta_MO(KT)>1
        Stab_T(KT)=5;
        Stab_m(KT)=5;
    else
        Stab_T(KT)=5*Zeta_MO(KT);
        Stab_m(KT)=5*Zeta_MO(KT);
    end
end

Velo_fric(KT)=U(KT)*VK_Const/(log((z_ref-d0_disp+z_srm)/z_srm)+Stab_m(KT));

Resis_a(KT)=(log((z_ref-d0_disp+z_srT)/z_srT)+Stab_T(KT))/(VK_Const*Velo_fric(KT));     %(s.cm^-1)

Resis_s(KT)=10*exp(35.63*(0.205-Theta_LL_sur(KT)))/100; %(-805+4140*(Theta_s(J)-Theta_LL_sur(KT)))/100;  % s.m^-1----->0.001s.cm^-1

r_a_SOIL(KT) = log((2.0)/0.001)*log(2.0/0.001)/((k^2)*U(KT))*100;   %(s.m^-1) ORIGINAL Zom=0.0058
Evapo(KT)=(RHOV_sur(KT)-RHOV_A(KT))/(Resis_s(KT)+r_a_SOIL(KT)/100);

% PT/PE - Penman-Montheith
% mm.day-1
% FAO56 pag19 eq3
% VEG
PT_PM_VEG(KT) = (DELTA(KT)*(Rn(KT))+3600*ro_a(KT)*Cp*(e0_Ta(KT)-e_a(KT))/r_a_VEG(KT))/(lambdav*(DELTA(KT) + gama*(1+r_s_VEG(KT)/r_a_VEG(KT))))/3600;

if LAI(KT)==0 || hh_v(KT)==0
    PT_PM_VEG(KT)=0;
end
PE_PM_SOIL(KT) = (DELTA(KT)*(Rn_SOIL(KT))+3600*ro_a(KT)*Cp*(e0_Ta(KT)-e_a(KT))/r_a_SOIL(KT))/(lambdav*(DELTA(KT) + gama*(1+r_s_SOIL(KT)/r_a_SOIL(KT))))/3600;
Evap(KT)=0.1*PE_PM_SOIL(KT); % transfer to second value-G_SOIL(KT)
EVAP(KT,1)=Evap(KT);
Tp_t(KT)=0.1*PT_PM_VEG(KT); %transfer to second value
TP_t(KT,1)=Tp_t(KT);

if rwuef==1
    if KT<=3288+1103
        H1=0;H2=-31;H4=-8000;H3L=-600;H3H=-300;
        if Tp_t(KT)<0.02/3600
            H3=H3L;
        elseif Tp_t(KT)>0.05/3600
            H3=H3H;
        else
            H3=H3H+(H3L-H3H)/(0.03/3600)*(0.05/3600-Tp_t(KT));
        end
    else
        H1=-1;H2=-5;H4=-16000;H3L=-600;H3H=-300;
        if Tp_t(KT)<0.02/3600
            H3=H3L;
        elseif Tp_t(KT)>0.05/3600
            H3=H3H;
        else
            H3=H3H+(H3L-H3H)/(0.03/3600)*(0.05/3600-Tp_t(KT));
        end
    end
    % piecewise linear reduction function
    MN=0;
    for ML=1:NL
        for ND=1:2
            MN=ML+ND-1;
            if hh(MN) >=H1
                alpha_h(ML,ND) = 0;
            elseif  hh(MN) >=H2
                alpha_h(ML,ND) = (H1-hh(MN))/(H1-H2);
            elseif  hh(MN) >=H3
                alpha_h(ML,ND) = 1;
            elseif  hh(MN) >=H4
                alpha_h(ML,ND) = (hh(MN)-H4)/(H3-H4);
            else
                alpha_h(ML,ND) = 0;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% Root lenth distribution %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Elmn_Lnth=0;
    RL=300;
    Elmn_Lnth(1)=0;
    RB=0.9;
    LR=log(0.01)/log(RB);
    RY=1-RB^(LR);
    
    if LR<=1
        for ML=1:NL-1      % ignore the surface root water uptake 1cm
            for ND=1:2
                MN=ML+ND-1;
                bx(ML,ND)=0;
            end
        end
    else
        for ML=2:NL
            Elmn_Lnth(ML)=Elmn_Lnth(ML-1)+DeltZ(ML-1);
            if Elmn_Lnth<RL-LR      %(KT)
                FRY(ML)=1;
            else
                FRY(ML)=(1-RB^(RL-Elmn_Lnth(ML)))/RY;
            end
        end
        for ML=1:NL-1
            bx(ML)=FRY(ML)-FRY(ML+1);
            if bx(ML)<0
                bx(ML)=0;
            end
            bx(NL)=FRY(NL);
        end
        for ML=1:NL
            for ND=1:2
                MN=ML+ND-1;
                bx(ML,ND)=bx(MN);
            end
        end
    end
    %root zone water uptake
    Trap_1(KT)=0;
    for ML=1:NL
        for ND=1:2
            MN=ML+ND-1;
            Srt_1(ML,ND)=alpha_h(ML,ND)*bx(ML,ND)*Tp_t(KT);
            if TT(ML)<0
                Srt_1(ML:NL,ND)=0;
            end
        end
        Trap_1(KT)=Trap_1(KT)+(Srt_1(ML,1)+Srt_1(ML,2))/2*DeltZ(ML);   % root water uptake integration by DeltZ;
    end
    
    %     % consideration of water compensation effect
    if Tp_t(KT)==0
        Trap(KT)=0;
    else
        wt(KT)=Trap_1(KT)/Tp_t(KT);
        wc=1; % compensation coefficient
        Trap(KT)=0;
        if wt(KT)<wc
            for ML=1:NL
                for ND=1:2
                    MN=ML+ND-1;
                    Srt(ML,ND)=alpha_h(ML,ND)*bx(ML,ND)*Tp_t(KT)/wc;
                    if TT(ML)<0
                        Srt(ML:NL,ND)=0;
                    end
                end
                Trap(KT)=Trap(KT)+(Srt(ML,1)+Srt(ML,2))/2*DeltZ(ML);   % root water uptake integration by DeltZ;
            end
        else
            for ML=1:NL
                for ND=1:2
                    MN=ML+ND-1;
                    Srt(ML,ND)=alpha_h(ML,ND)*bx(ML,ND)*Tp_t(KT)/wt(KT);
                end
                Trap(KT)=Trap(KT)+(Srt(ML,1)+Srt(ML,2))/2*DeltZ(ML);   % root water uptake integration by DeltZ;
            end
        end
    end
end