% Calculates Energy and Mass Flux at any instant

function [FM,Q,QM,MR,QE,E,TSURF,QH,QNET]=QFM(UB,W,A,TA,PRAIN,PS,WS,RH,QSI,QLI,RKN,IRADFL,...
                     QNETOB,CDH,CDE,RID,param,sitev)
                     
%  Site variables
      FC=sitev(1);     % Forest cover fraction (0-1)
%     df=sitev(2)     % Drift factor
      PR=sitev(3);     % Atmospheri%Pressure (Pa)
      QG=sitev(4);     % Ground heat flux (KJ/m^2/hr)  This is more logically an
%     input variable, but is put here because it is never known at each
%     time step.  Usually it will be assigned a value 0.

% Parameters
%     TR=param(1);     % Temperature above which all is rain (3 C)
%     TS=param(2);     % Temperature below which all is snow (-1 C)
      TO=param(3);     % Temperature of freezing (0 C)
      TK=param(4);     % Temperature to convert C to K (273.15)
      ES=param(5);     % emmissivity of snow (nominally 0.99)
      SBC=param(6);    % Stefan boltzman constant (2.0747e-7 KJ/m^2/hr/K)
      HF =param(7);    % Heat of fusion (333.5 KJ/kg)
      HNEU=param(8);   % Heat of Vaporization (Ice to Vapor, 2834 KJ/kg)
      CW =param(9);    % Water Heat Capacity (4.18 KJ/kg/C)
      CS =param(10);   % Ice heat capacity (2.09 KJ/kg/C)
      CG =param(11);   % Ground heat capacity (nominally 2.09 KJ/kg/C)
      CP=param(12);    % Air Heat Capacity (1.005) (KJ/kg/K)
      RA =param(13);   % Ideal Gas constant for dry air (287 J/kg/K)
      K=param(14);     % Von Karmans constant (0.4)
      Z=param(15);     % Nominal measurement height for air temperature and humidity (2m)
      Zo=param(16);    % Surface aerodynamic roughness (m)
      hff=param(17);   % Factor to convert /s into /hr (3600)
      RHOI=param(18);  % Density of Ice (917 kg/m^3)
      RHOW=param(19);  % Density of Water (1000 kg/m^3)
      RHO=param(20);   % Snow Density (Nominally 450 kg/m^3)
      RHOG=param(21);  % Soil Density (nominally 1700 kg/m^3)
      LC=param(22);    % Liquid holding capacity of snow (0.05)
      KS=param(23);    % Snow Saturated hydraulic conductivity (160 m/hr)
      DE=param(24);    % Thermally active depth of soil (0.4 m)
      RS=param(25);    % Snow Surface thermal conductance (m/hr)
      G=param(26);     % Gravitational acceleration (9.81 m/s^2)
      fstab=param(30);  % Stability correction control parameter
      EA=SVPW(TA)*RH;   %The saturation vapour pressure over water is used 
%      because most instruments report relative humidity relative to water.
       [QP]=QPF(PRAIN,TA,TO,PS,RHOW,HF,CW,CS);
       [TAVE]=TAVG(UB,W,RHOW,CS,TO,RHOG,DE,CG,HF);
       [TSURF] = SRFTMP(QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,...
                 RKN,HNEU,ES,SBC,CS,RS,W,QNETOB,IRADFL,WS,Z,G,FC,...
                 fstab);
       QLE=(1-FC)*ES*SBC*(TSURF+TK)^4; % outgoing longwave radiation Eq. 26
       QLNET=QLI-QLE;                  % Net longwave radiation
       [QH,QE,E]=TURBFLUX(PR,RA,TA,TK,TSURF,Z,G,CP,RKN,WS,EA,...
                         RHOW,HNEU,fstab);
       [MR]=FMELT(UB,RHOW,W,HF,LC,RID,KS,PRAIN);
       if MR<=0
           MR=0;
       end
%       MR in m/hr
       QM=MR*RHOW*(HF+(TAVE-TO)*CW);        % Include advection of
%        meltwater/rain that is at tave so that the model does not
%        crash when there is no snow and it rains.
%       QM in kj/m2/hr 
       if (IRADFL==0) 
          QNET = QSI*(1.0-A)+QLNET;
       else
          QNET = QNETOB;
       end
       Q = QNET+QP+QG+QH+QE-QM;       
       FM=PRAIN+PS-MR-E;
       
end
