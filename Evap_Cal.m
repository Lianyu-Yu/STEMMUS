global RHOV_sur RHOV_A Resis_a Resis_s P_Va Velo_fric Theta_LL_sur  % RHOV_sur and Theta_L_sur should be stored at each time step.
global z_ref z_srT z_srm VK_Const d0_disp U_wind MO Ta U Ts Zeta_MO Stab_m Stab_T       % U_wind is the mean wind speed at height z_ref (m��s^-1), U is the wind speed at each time step.
global Rv g HR_a NL NN Evap KT RHOV Theta_LL EVAP
global Evaptranp_Cal DayNum t Ep DURTN  ETp Tp Tao LAI Tp_t Trap AFTP_TIME
global H1 H2 H3 H4 alpha_h bx LR Lm fr RL0 Srt Elmn_Lnth DeltZ RL TIME rwuef hh
global PT_PM_0 PT_PM_VEG PE_PM_SOIL T_act Theta_s J Resis_s1 Resis_s2 Resis_s3 Resis_s4 Resis_s5 Resis_a1
global Rns Rns_SOIL RnL G_SOIL r_s_VEG r_s_SOIL r_a_VEG r_a_SOIL Rn_SOIL Rn  Rs Ra w w1 w2 ws Rs0 e0_Ts e_a_Ts e0_Ta e_a wt wc Srt_1 rl LAI_act rl_min
global Es PME Kcb ET Theta_r Coefficient_n Coefficient_Alpha

%%%%%%% LAI and light extinction coefficient calculation %%%%%%%%%%%%%%%%%%
AFTP_TIME=TIME/86400+9;   %9 is the daynumber initial;

if AFTP_TIME<14
    LAI(KT)=0; % emergance daynumber is 8
elseif AFTP_TIME<22
    LAI(KT)=(AFTP_TIME-14)*0.45/8; % emergance daynumber is 8
else
    LAI(KT)=-0.0021*AFTP_TIME^2+0.299*AFTP_TIME-5.1074;
end
if LAI(KT)<0
    LAI(KT)=0;
end
LAI(KT)=LAI(KT);
if LAI(KT)<=2
    LAI_act(KT)=LAI(KT);
elseif LAI(KT)<=4
    LAI_act(KT)=2;
else
    LAI_act(KT)=0.5*LAI(KT);
end
% LAI_act(KT)=0.67*LAI(KT);
Tao=0.56;  %light attenual coefficient

if Evaptranp_Cal==1

    % input data
    n(J)=Coefficient_n(J);
    Alpha(J)=Coefficient_Alpha(J);
    m(J)=1-1/n(J);
    AFTP_TIME=TIME/86400+27;   %27 is the daynumber initial;
    Theta_LL_sur1(KT)=Theta_LL(NL,2);
    Theta_LL_sur2(KT)=Theta_LL(NL-14,2);
    Theta_LL_sat(KT)=Theta_r(J)+(Theta_s(J)-Theta_r(J))/(1+abs(Alpha(J)*200)^n(J))^m(J);
    coef_e=0.9; % 0.8-1.0 Johns 1982, Kemp 1997
    coef_p=2.15;  %2-2.3
    Kcbmax=1.20;  %for maize 1.10, for wheat 1.07 (allen 2009;duchemin 2006)
    coef_kd=-0.7; %-0.84 for wheat
    Kcb(KT)=Kcbmax*(1-exp(coef_kd*LAI(KT)));
    if TIME<DURTN
        DayNum=fix(TIME/3600/24)+1;
        t=TIME-(DayNum-1)*86400;
        ETp=0.1.*PME(14:115);
        ET=0.1.*PT_PM_0(14:115);
        Ep(KT)=(exp(-1*(Tao*LAI(KT))))*ETp(DayNum);
        
        if hh(NN)>-1e5
            if (Theta_LL_sur1(KT)/Theta_LL_sat(KT))>((Ep(KT)/coef_e)^0.5)
                Es(KT)=Ep(KT);
            else
                Es(KT)=coef_e*(Theta_LL_sur1(KT)/Theta_LL_sat(KT))^coef_p;
            end
        else
            Es(KT)=coef_e*((Theta_LL_sur1(KT)+Theta_LL_sur2(KT))/2/Theta_LL_sat(KT))^coef_p;
        end
        % generate E and T function with time
        if t>0.264*24*3600 && t<0.736*24*3600
            
            Tp(KT)=Kcb(KT)*ET(DayNum); % Tao LAI set as constant
            Evap(KT)=(2.75*sin(2*pi()*t/3600/24-pi()/2)/86400)*Es(KT); % transfer to second value
            EVAP(KT,1)=Evap(KT);
            Tp_t(KT)=(2.75*sin(2*pi()*t/3600/24-pi()/2)/86400)*Tp(KT); % transfer to second value
            TP_t(KT,1)=Tp_t(KT);
        else
            Tp(KT)=Kcb(KT)*ET(DayNum); % Tao LAI set as constant
            %           Tp(KT)=(1-exp(-1*(Tao*LAI(KT))))*ET(DayNum); % Tao LAI set as constant
            Evap(KT)=(0.24/86400)*Es(KT); % transfer to second value
            EVAP(KT,1)=Evap(KT);
            Tp_t(KT)=(0.24/86400)*Tp(KT); % transfer to second value
            TP_t(KT,1)=Tp_t(KT);
        end
    else
        DayNum=fix(TIME/3600/24);
        t=TIME-(DayNum-1)*86400;
        ETp=0.1.*PME(14:115);
        ET=0.1.*PT_PM_0(14:115);
        Ep(KT)=(exp(-1*(Tao*LAI(KT))))*ETp(DayNum);
        %        Kcb(KT)=Kcbmax*(1-exp(coef_kd*LAI_act(KT)));
        if hh(NN)>-1e5
            if (Theta_LL_sur1(KT)/Theta_LL_sat(KT))>((Ep(KT)/coef_e)^0.5)
                Es(KT)=Ep(KT);
            else
                Es(KT)=coef_e*(Theta_LL_sur1(KT)/Theta_LL_sat(KT))^coef_p;
            end
        else
            Es(KT)=coef_e*((Theta_LL_sur1(KT)+Theta_LL_sur2(KT))/2/Theta_LL_sat(KT))^coef_p;
        end
        % generate E and T function with time
        if t>0.264*24*3600 && t<0.736*24*3600
            
            Tp(KT)=Kcb(KT)*ET(DayNum); % Tao LAI set as constant
            Evap(KT)=(2.75*sin(2*pi()*t/3600/24-pi()/2)/86400)*Es(KT); % transfer to second value
            EVAP(KT,1)=Evap(KT);
            Tp_t(KT)=(2.75*sin(2*pi()*t/3600/24-pi()/2)/86400)*Tp(KT); % transfer to second value
            TP_t(KT,1)=Tp_t(KT);
        else
            Tp(KT)=Kcb(KT)*ET(DayNum); % Tao LAI set as constant
            %           Tp(KT)=(1-exp(-1*(Tao*LAI(KT))))*ET(DayNum); % Tao LAI set as constant
            Evap(KT)=(0.24/86400)*Es(KT); % transfer to second value
            EVAP(KT,1)=Evap(KT);
            Tp_t(KT)=(0.24/86400)*Tp(KT); % transfer to second value
            TP_t(KT,1)=Tp_t(KT);
        end
    end
else
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
    Z=521;             % altitute of the location(m)
    as=0.25;           % regression constant, expressKTg the fraction of extraterrestrial radiation FAO56 pag50
    bs=0.5;
    alfa=0.23;         % albeo of vegetation set as 0.23
    z_m=10;            % observation height of wKTd speed; 10m
    Lz=240*pi()/180;   % latitude of Beijing time zone west of Greenwich
    Lm=252*pi()/180;    % latitude of Local time, west of Greenwich
    % albedo of soil calculation;
    Theta_LL_sur(KT)=Theta_LL(NL,2);
    if Theta_LL_sur(KT)<0.1
        alfa_s(KT)=0.25;
    elseif Theta_LL_sur(KT)<0.25
        alfa_s(KT)=0.35-Theta_LL_sur(KT);
    else
        alfa_s(KT)=0.1;
    end
    JN(KT)=fix(TIME/3600/24)+174;    % day number
    n=[4.8	0	10.6	2.9	13	12.3	9.3	10.9	2.6	0	4.1	0	11.9	11.4	8.1	0	0	0	4.7	5.4	10.2	0	0	12.1	0	0	5.4	11.1	0	1.5	0.7	0	2.5	8.7	6	3.9	0	1.2	10.5	8.6	7.3	8.8	9.8	10.8	8.6	0	4.6	8.9	3.2	3	9	7.9	4	7.2	6.3	5.1	9.2	8.9	9.7	8.2	4.5	3.1	0	4.9	8.1	0	0	11.6	11.2	7.7	7.2	0	2.9	0	3	9	0	0	8	9.1	4.5	4.7	11	11.2	9.7	8.8	7.3	0	0	0	4.4	0	0	3.6	0	0	8.6	0.8	8.6	0	7.4	3.1];
    h_v=[0.047333333	0.071	0.094666667	0.118333333	0.142	0.165666667	0.189333333	0.213	0.236666667	0.260333333	0.284	0.307666667	0.331333333	0.355	0.390909091	0.426818182	0.462727273	0.498636364	0.534545455	0.570454545	0.606363636	0.642272727	0.678181818	0.714090909	0.75	0.807142857	0.864285714	0.921428571	0.978571429	1.035714286	1.092857143	1.15	1.21625	1.2825	1.34875	1.415	1.48125	1.5475	1.61375	1.68	1.715	1.75	1.785	1.82	1.855	1.89	1.925	1.96	1.995	2.03	2.065	2.1	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135	2.135];
    rl_min=[139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	121	121	121	121	121	121	121	121	121	121	121	121	121	121	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	139	239	239	239	239	239	239];
    DayNum=fix(TIME/3600/24)+1;
    n(KT)=n(DayNum);
    h_v(KT)=h_v(DayNum);
    %rl_min(KT)=139;
    
    rl_min(KT)=rl_min(DayNum);
    % 6-23 TO 7-31
    %Kcb=Mdata(:,15);
    
    % Calculation procedure
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
    % compute w1 - solar time angle at the beginning of the period (time)
    % [rad]
    % FAO56 pag47 Eq29
    tl = 1;  %  hourly data
    w1(KT) = (w(KT) - pi()*tl/24.0);
    
    % compute w2 - solar time angle at the end of the period (time + 1h)
    % [rad]
    % FAO56 pag47 Eq30
    w2(KT) = w(KT) + pi()*tl/24.0;
    
    % compute ws - sunset hour angle
    % [rad]
    % FAO56 pag47 Eq25
    ws(KT)=acos((-1)*tan(0.599)*tan(delta(KT)));  %for daily duration
    
    % compute Ra - extraterrestrial radiation
    % [MJ.m-2.hour-1]
    % FAO56 pag47 Eq28
    if w(KT)> -ws(KT) &&  w(KT) < ws(KT)
        Ra(KT)=12*60/pi()*Gsc*dr(KT)*((w2(KT)-w1(KT))*sin(0.599)*sin(delta(KT)) + cos(0.599)*cos(delta(KT))*(sin(w2(KT))-sin(w1(KT))));
    else
        Ra(KT)=0;
    end
    % compute Rs0 - clear-sky solar (shortwave) radiation
    % [MJ.m-2.hour-1]
    % FAO56 pag51 Eq37
    Rs0(KT) = (0.75+2E-5*Z)*Ra(KT);
    %    Rs0_Watts = Rs0*24.0/0.08864
    % daylight hours N
    N(KT)=24*ws(KT)/pi();
    
    % compute Rs - SHORTWAVE RADIATION
    % [MJ.m-2.hour-1]
    % FAO56 pag51 Eq37
    Rs(KT)=(as+bs*n(KT)/N(KT))*Ra(KT);
    
    % compute Rns - NET SHORTWAVE RADIATION
    % [MJ.m-2.day-1]
    % FAO56 pag51 Eq37
    % for each type of vegetation, crop and soil (albedo dependent)
    Rns(KT)= (1-alfa)*Rs(KT);
    Rns_SOIL(KT) = (1 - alfa_s(KT))*Rs(KT);
    % compute Rnl - NET LONGWAVE RADIATION
    % [MJ.m-2.hour-1]
    % FAO56 pag51 Eq37 and pag74 of hourly computKTg
    r_sunset=[];
    r_angle=[];
    R_i=[];
    if (ws(KT) - 0.52) <= w(KT) && w(KT)<= (ws(KT) - 0.10) %FAO56: (ws(KT) - 0.79) <= w(KT) <= (ws(KT) - 0.52)
        R_i = 1;
        if Rs0(KT) > 0
            if Rs(KT)/Rs0(KT) > 0.3
                r_sunset = Rs(KT)/Rs0(KT);
            else
                r_sunset = 0.3;
            end
        else
            r_sunset = 0.75;  % see FAO56 pag75
        end
    end
    if (ws(KT) - 0.10) < w(KT) || w(KT) <= (-ws(KT)+ 0.10)
        if R_i>0
            r_angle(KT)=r_sunset;
        else
            r_angle(KT)=0.75;  %see FAO56 pag75
        end
    else
        r_angle(KT)=Rs(KT)/Rs0(KT);
    end
    RnL(KT)=(sigma/24*((Ta(KT) + 273.16)^4)*(0.34-0.14*sqrt(e_a(KT)))*(1.35*r_angle(KT)-0.35));
    
    if RnL(KT)<0
        r_angle(KT)=0.8;
        RnL(KT)=(sigma/24*((Ta(KT) + 273.16)^4)*(0.34-0.14*sqrt(e_a(KT)))*(1.35*r_angle(KT)-0.35));
    end
    
    Rn(KT) = Rns(KT) - RnL(KT);  % net radiation for vegetation
    %    Rn_SOIL(KT) = Rns_SOIL(KT) - RnL(KT);  % net radiation for vegetation
    Rn_SOIL(KT) =Rn(KT)*exp(-1*(Tao*LAI(KT)));  % net radiation for soil
    % soil heat flux
    t=TIME-(fix(TIME/3600/24))*86400;
    if t>0.264*24*3600 && t<0.736*24*3600
        G(KT)=0.1*Rn(KT);
        G_SOIL(KT)=0.1*Rn_SOIL(KT);
    else
        G(KT)=0.5*Rn(KT);
        G_SOIL(KT)=0.5*Rn_SOIL(KT);
    end
    
    %% SURFACE RESISTANCE PARAMETERS CALCULATION
    R_a=0.81;R_b=0.004*24*11.6;R_c=0.05;
    % R_fun(KT)=((R_b*Rns(KT)+R_c)/(R_a*(R_b*Rns(KT)+1)));
    rl(KT)=rl_min(KT)/((R_b*Rns(KT)+R_c)/(R_a*(R_b*Rns(KT)+1)));
    
    % r_s - SURFACE RESISTANCE
    % [s.m-1]
    % VEG: Dingman pag 208 (canopy conductance) (equivalent to FAO56 pag21 Eq5)
    r_s_VEG(KT) = rl(KT)/LAI_act(KT);
    
    % SOIL: equation 20 of van de Griend and Owe, 1994
    %Theta_LL_sur(KT)=Theta_LL(NL,2);
    
    r_s_SOIL(KT)=10.0*exp(0.3563*100.0*(0.25-Theta_LL_sur(KT)));   % 0.25 set as minmum soil moisture for potential evaporation
    %r_s_SOIL(i)=10.0*exp(0.3563*100.0*(fc(i)-por(i)));
    % correction wKTdspeed measurement and scalKTg at h+2m
    % [m.s-1]
    % FAO56 pag56 eq47
    
    % r_a - AERODYNAMIC RESISTANCE
    % [s.m-1]
    % FAO56 pag20 eq4- (d - zero displacement plane, z_0m - roughness length momentum transfer, z_0h - roughness length heat and vapour transfer, [m], FAO56 pag21 BOX4
    r_a_VEG(KT) = log((2-(2*h_v(KT)/3))/(0.123*h_v(KT)))*log((2-(2*h_v(KT)/3))/(0.0123*h_v(KT)))/((k^2)*U(KT))*100;
    % r_a of SOIL
    % Liu www.hydrol-earth-syst-sci.net/11/769/2007/
    % equation for neutral conditions (eq. 9)
    % only function of ws, it is assumed that roughness are the same for any type of soil
    RHOV_sur(KT)=RHOV(NN);
    Theta_LL_sur(KT)=Theta_LL(NL,2);
    
    P_Va(KT)=0.611*exp(17.27*Ta(KT)/(Ta(KT)+237.3))*HR_a(KT);  %The aTaospheric vapor pressure (KPa)  (1000Pa=1000.1000.g.100^-1.cm^-1.s^-2)
    
    RHOV_A(KT)=P_Va(KT)*1e4/(Rv*(Ta(KT)+273.15));              %  g.cm^-3;  Rv-cm^2.s^-2.Cels^-1
    
    z_ref=200;          % cm The reference height of tempearature measurement (usually 2 m)
    d0_disp=0;          % cm The zero-plane displacement (=0 m)
    z_srT=0.1;          % cm The surface roughness for the heat flux (=0.001m)
    VK_Const=0.41;   % The von Karman constant (=0.41)
    z_srm=0.1;          % cm The surface roughness for momentum flux (=0.001m)
    U_wind=198.4597; % The mean wKTd speed at reference height.(cm.s^-1)
    
    MO(KT)=((Ta(KT)+273.15)*U(KT)^2)/(g*(Ta(KT)-Ts(KT))*log(z_ref/z_srm));  % Wind speed should be KT cm.s^-1, MO-cm;
    
    Zeta_MO(KT)=z_ref/MO(KT);
    
    if abs(Ta(KT)-Ts(KT))<=0.01
        Stab_m(KT)=0;
        Stab_T(KT)=0;
    elseif Ta(KT)<Ts(KT) || Zeta_MO(KT)<0
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
    
    Resis_a(KT)=((log((z_ref-d0_disp+z_srT)/z_srT)+Stab_T(KT))/(VK_Const*Velo_fric(KT)))*100;     %(s.cm^-1)
    r_a_SOIL(KT) = log((2.0)/0.0058)*log(2.0/0.0058)/((k^2)*U(KT))*100;   %(s.m^-1)
    
    % PT/PE - Penman-Montheith
    % mm.day-1
    % FAO56 pag19 eq3
    % VEG
    PT_PM_VEG(KT) = (DELTA(KT)*(Rn(KT))+3600*ro_a(KT)*Cp*(e0_Ta(KT)-e_a(KT))/r_a_VEG(KT))/(lambdav*(DELTA(KT) + gama*(1+r_s_VEG(KT)/r_a_VEG(KT))))/3600;
    % reference et ET0
    %PT_PM_0(KT) = (0.408*DELTA(KT)*Rn(KT)+gama*900/(Ta(KT)+273)*(e0_Ta(KT)-e_a(KT))*u_2(KT))/(DELTA(KT) + gama*(1+0.34*u_2(KT)));
    %T_act(KT)=PT_PM_0(KT)*Kcb(KT);
    % for SOIL
    PE_PM_SOIL(KT) = (DELTA(KT)*(Rn_SOIL(KT)-G_SOIL(KT))+3600*ro_a(KT)*Cp*(e0_Ta(KT)-e_a(KT))/r_a_SOIL(KT))/(lambdav*(DELTA(KT) + gama*(1+r_s_SOIL(KT)/r_a_SOIL(KT))))/3600;
    Evap(KT)=0.1*PE_PM_SOIL(KT); % transfer to second value
    EVAP(KT,1)=Evap(KT);
    Tp_t(KT)=0.1*PT_PM_VEG(KT); % transfer to second value
    TP_t(KT,1)=Tp_t(KT);
    
end
if rwuef==1
    % water stress function parameters
    H1=-15;H2=-50;H4=-9000;H3L=-900;H3H=-500;
    if Tp_t(KT)<0.02/3600
        H3=H3L;
    elseif Tp_t(KT)>0.05/3600
        H3=H3H;
    else
        H3=H3H+(H3L-H3H)/(0.03/3600)*(0.05/3600-Tp_t(KT));
    end
    
    % piecewise linear reduction function
    MN=0;
    for ML=1:NL
        for ND=1:2
            MN=ML+ND-1;
            if hh(MN) >=H1,
                alpha_h(ML,ND) = 0;
            elseif  hh(MN) >=H2,
                alpha_h(ML,ND) = (H1-hh(MN))/(H1-H2);
            elseif  hh(MN) >=H3,
                alpha_h(ML,ND) = 1;
            elseif  hh(MN) >=H4,
                alpha_h(ML,ND) = (hh(MN)-H4)/(H3-H4);
            else
                alpha_h(ML,ND) = 0;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% Root lenth distribution %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Lm=120;
    RL0=20;   
    r=9.48915E-07;  % root growth rate cm/s
    fr(KT)=RL0/(RL0+(Lm-RL0)*exp((-1)*(r*TIME)));
    LR(KT)=Lm*fr(KT);
    RL=300;
    Elmn_Lnth=0;
    
    if LR(KT)<=1
        for ML=1:NL-1      % ignore the surface root water uptake 1cm
            for ND=1:2
                MN=ML+ND-1;
                bx(ML,ND)=0;
            end
        end
    else
        for ML=1:NL
            Elmn_Lnth=Elmn_Lnth+DeltZ(ML);
            if Elmn_Lnth<RL-LR(KT)
                bx(ML)=0;
            elseif Elmn_Lnth>=RL-LR(KT) && Elmn_Lnth<RL-0.2*LR(KT)
                bx(ML)=2.0833*(1-(RL-Elmn_Lnth)/LR(KT))/LR(KT);
            else
                bx(ML)=1.66667/LR(KT);
            end
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