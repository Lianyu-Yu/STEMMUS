function [statevu,statevw,cump,cume,cummr,outv]=snow_Calc(Delt_t,Ta,U,HR_a,Rns,Rnl,KT,Rn,Precip)
global ALBD_G coszen1 coszen0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
io=4914.0;    %solar constant  Kj/m^2/hr
param=[3.5 0 0 273.15 0.99 2.0747e-7...
    333.5 2834 4.18 2.09 2.09 1.005...
    287 0.4 2 0.005 3600 ...
    917 1000 450 1700 ...
    0.05 20 0.4 0.02 9.81 0.25 0.85 0.65 0];

sitev=[0.0  1 67641.648   0 0.0001 ...     %fc  df   pr(PA)  qg  aep
    0 0  34];  %slope aspect lat
statev=[0 0.0000 0.267];%A(1,6:8);
SLOPE=sitev(6);
AZI=sitev(7);
LAT=sitev(8);
% set control flags
iflag(1)=2;%IRADFL;   %radiation
iflag(2)=1;        %Printing
iflag(3)=9;        %Output unit
iflag(4)=1;        %Albedo calculation

%  initialize variables for mass balance
W1=statev(2);
cump=0;
cume=0;
cummr=0;
fac=0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsv=5;npar=30;nxv=3;niv=7;nov=13;
DT=Delt_t/3600;   % time step (hr)
YEAR=2015;
MONTH=12;
DAY=1;
HOUR=9.5;
IRAD=2;
% State variables - These serve as initial conditions
UB=78.5;%vub(KT-1);%statev(1);    % Snow Energy Content  (KJ/m^2)
W=0.000001;%vw(KT-1);%statev(2);     % Snow Water Equivalent (m) relative to T = 0 C solid phase

%%%%%%%%%%%%%%%   Site variables    %%%%%%%%%%%%%%%%%%%%
FC=sitev(1);     %  Forest cover fraction (0-1)
DF=sitev(2);     %  Drift factor
PR=sitev(3);     %  Atmospheric Pressure (Pa)
QG=sitev(4);     %  Ground heat flux (KJ/m^2/hr)  This is more logically an
%             input variable,but is put here because it is never known at each
%                 time step. Usually it will be assigned a value 0.
aep=sitev(5);   %  Albedo extinction parameter to smooth
%    transition of albedo when snow is shallow. Depends on Veg. height (m)
SLOPE=sitev(6);
AZI=sitev(7);
LAT=sitev(8);
%%%%%%%%%%%%%%%%%%%    Parameters    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TR=3.5;%param(1);     %Temperature above which all is rain (3 C)
TS=0;%param(2);     %Temperature below which all is snow (-1 C)
TO=param(3);     %Temperature of freezing (0 C)
TK=param(4);     %Temperature to convert C to K (273.15)
ES=param(5);     %emmissivity of snow (nominally 0.99)
SBC=param(6);    %Stefan boltzman constant (2.0747e-7 KJ/m^2/hr/K)
HF =param(7);    %Heat of fusion (333.5 KJ/kg)
HNEU=param(8);   %Heat of Vaporization (Ice to Vapor, 2834 KJ/kg)
CW =param(9);    %Water Heat Capacity (4.18 KJ/kg/C)
CS =param(10);   %Ice heat capacity (2.09 KJ/kg/C)
CG =param(11);   %Ground heat capacity (nominally 2.09 KJ/kg/C)
CP=param(12);    %Air Heat Capacity (1.005 KJ/kg/K)
RA =param(13);   %Ideal Gas constant for dry air (287 J/kg/K)
k=param(14);     %Von Karmans constant (0.4)
Z=param(15);     %Nominal meas. height for air temp. and humidity (2m)
Zo=param(16);    %Surface aerodynamic roughness (m)0.001;%
hff=param(17);   %Factor to convert /s into /hr (3600)
RHOI=param(18);  %Density of Ice (917 kg/m^3)
RHOW=param(19);  %Density of Water (1000 kg/m^3)
RHO=param(20);   %Snow Density (Nominally 450 kg/m^3)
RHOG=param(21);  %Soil Density (nominally 1700 kg/m^3)
LC=param(22);    %Liquid holding capacity of snow (0.05)
KS=param(23);    %Snow Saturated hydraulic conductivity (160 m/hr)
DE=param(24);    %Thermally active depth of soil (0.4 m)
RS=param(25);    %Snow Surface thermal conductance (m/hr)
G=param(26);     %Gravitational acceleration (9.81 m/s^2)
abg=ALBD_G(KT);  %param(27);   %Bare ground albedo  (0.25)
avo=param(28);   %Visual new snow albedo (0.95)
anir0 = param(29); % NIR new snow albedo (0.65)
fstab = param(30); % Stability correction control parameter
% 0 = no corrections, 1 = full corrections
bca=0.8;    % values Bristow and Campbell
bcc=2.4;    % values Bristow and Campbell
dtbar(1:12)=[5.712903 4.35 6.890322...
    8.660001 8.93871 10.01 9.541936 9.038710 7.160001 8.106450 5.92332 5.058064];

%%%%%%%%%%%%%%%%%   control flags  %%%%%%%%%%%%%%%%%%%%%
IRADFL=iflag(1);
pflag=iflag(2);
ounit=iflag(3);
% Calculate constants
%  These need only be calculated once.
%  The model is more efficient for large nt since it saves on these calc's
CDH=CP*(k/log(Z/Zo))^2*hff;   %  Factor in heat transfer
CDE=0.622*HNEU*(k/log(Z/Zo))^2*hff;  % Factor in vapor transfer
CD=k*k*hff/(log(Z/Zo)^2)*(1-0.8*FC);   % factor in turbulent fluxes
%  The FC correction is applied here to hit sensible and latent heat fluxes
%  and evaporation rate in one place and maintain consistency in the surface
%   temperature and stability calculations. It is consistent with wind speed
%   reduction by a factor of 1-0.8 FC which is the most reasonable physical
%   justification for this approach.
%   I recognize that this is not a very good way to parameterize vegetation.
%   This will be an area of future improvements (I hope).
%    FC is also used below to adjust radiation inputs.
RRHOI=RHOI/RHOW;
RRHO=RHO/RHOW;
RID=1.0/RRHO-1.0/RRHOI;
%   Loop for each time step
%     for i = 1:nt
%   Input variables
TA=Ta(KT);%(1,i);    % Air temperature input (Degrees C)
P=Precip(KT)*36;     % Precipitation rate input (m/hr)
WS=U(KT)/100;    % Wind Speed (m/s)
RH=HR_a(KT);    % Relative humidity (fraction 0-1)
trange=0;    % temperature range
unit=1000;               % Rns in MJ/M2/hr
%         unit=8.64/24/100;          % Rns in W/M2
QSIOBS=Rns(KT)*unit;    % Observed short wave radiation UNIT convert TO KJ/M2/hr
QNETOBS=Rn(KT)*unit;    % Observed short wave radiation
QLIOBS=Rnl(KT)*unit;
[HRI,COSZEN]=hyri(MONTH,DAY,HOUR,DT,SLOPE,AZI,LAT);

if (IRAD<=1)
    [atff]=atf(trange,MONTH,dtbar,bca,bcc);
    if (IRAD==0)   %NEED TO ESTIMATE RADIATION BASED ON AIR TEMP.
        QSI=atff*io*HRI;  %INPT(5,1)
    else
        %     Need to call HYRI for horizontal surface to perform horizontal
        %     measurement adjustment
        [HRI0,COSZEN0]=hyri(MONTH,DAY,HOUR,DT,0,AZI,LAT);
        coszen0(KT)=COSZEN0;
        %            call hyri(YEAR,MONTH,DAY,HOUR,DT,0.0,AZI,LAT,HRI0,COSZEN)
        %     If HRI0 is 0 the sun should have set so QSIOBS should be 0.  If it is
        %     not it indicates a potential measurement problem. i.e. moonshine
        if (HRI0>0.0)
            QSI=QSIOBS*HRI/HRI0;
        else
            QSI=QSIOBS;
            if(QSIOBS ~= 0)
                print 'Warning you have nonzero nightime';
            end
        end
    end
    
    %       cloud cover fraction dgt 10/13/94
    cf = 1-atff/bca;
    [QLIFF]=qlif(TA,RH,TK,SBC,cf);
    %     call qlif(INPT(6,1),TA,RH,param(4),param(6),cf)
    QLI=QLIFF;
    IRADFL=0;
    QNETOB=0;
else
    IRADFL=1;
    QNETOB=QNETOBS;
    QSI=QSIOBS;
    QLI=QLIOBS;
end

%  UPDATETIME(YEAR,MONTH,DAY,HOUR,DT);
[YEAR,MONTH,DAY,HOUR]=UPDATETIME(YEAR,MONTH,DAY,HOUR,DT);

%  Separate rain and snow
PS=PARTSNOW(P,TA,TR,TS);
PRAIN = P - PS;
%  Increase precipitation as snow by drift multiplication factor
PS = PS*DF;
if (iflag(4)==1)
    %  Calculate albedo
    A=ALBEDO(statev(3),COSZEN,W/RRHO,aep,abg,avo,anir0);
    % Use of this albedo throughout time step neglects the
    %  changes due to new snow within a time step.
else
    A=statev(3);
end
%   Calculate neutral mass transfer coefficient
RKN=CD*WS;  % Eq. 36 P14
%   Adjust radiation inputs for the effect of forest canopy.
QSI=QSI*(1-FC);
QLI=QLI*(1-FC);
QNETOB=QNETOB*(1-FC);
%     FC corrections are also in the following subroutines.
%      QFM where outgoing longwave radiation is calculated.
%      SNOTMP where outgoing longwave radiation is used in the surface
%      temperature equilibrium approx.

[QH,QE,E,MR,QM,Q,FM,TSURF,QNET,W,UB]=PREDICORR(DT,UB,W,A,TA,PRAIN,PS,WS,RH,QSI,...
    QLI,IRADFL,RKN,QNETOB,CDH,CDE,RID,param,sitev);
%   Update snow surface age based on snowfall in time step
if(iflag(4)==1) %call agesn(statev(3),dt,ps,tsurf,tk)
    [statev(3)]=agesn(statev(3),DT,PS,TSURF,TK);
    %    accumulate for mass balance
    cump=cump+(PS+PRAIN)*DT;
    cume=cume+E*DT;
    cummr=cummr+MR*DT;
    tave=TAVG(UB,W,RHOW,CS,TO,RHOG,DE,CG,HF);   %  this call only
    %  Calculate albedo
    A=ALBEDO(statev(3),COSZEN,W/RRHO,aep,abg,avo,anir0);
end
coszen1(KT)=COSZEN;
statevu(KT,1)=UB;
statevw(KT,1)=W;
outv(KT,1)=PRAIN;    % Rainfall(m/h) Part of precipitation modeles as rain
outv(KT,2)=PS;       % Snowfall(m/h) Part of precipitation modeles as snow
outv(KT,3)=A;        % Albedo
outv(KT,4)=QH;       % Sensible heat flux Qh(KJ/m2/h)
outv(KT,5)=QE;       % Latent heat flux Qe(KJ/m2/h)
outv(KT,6)=E;        % Sublimation (m/h)
outv(KT,7)=MR;       % Melt outflow rate (m/h)
outv(KT,8)=QM;       % Heat advected by melt outflow(KJ/m2/h)
outv(KT,9)=Q;        % Total surface energy flux into the snow(KJ/m2/h)
outv(KT,10)=FM;      % Combined mass fluxes(dW/dt)(m/h)
outv(KT,11)=tave;    % Snow and underlying soil average temperature T (oC)
outv(KT,12)=TSURF;   % Snow surface temperature, Ts (oC)
outv(KT,13)=QNET;    % Net radiation (KJ/m2/h)
outv(KT,14)=QSI;      % Net shortwave radiation (KJ/m2/h)
outv(KT,15)=QLI;      % Net longwave radiation (KJ/m2/h)
%     end
end