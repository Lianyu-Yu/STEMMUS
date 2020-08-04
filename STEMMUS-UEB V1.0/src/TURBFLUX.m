% Calculates the turbulent heat fluxes (sensible and latent
%   heat fluxes) and condensation/sublimation.

function [QH,QE,E]=TURBFLUX(PR,RA,TA,TK,TS,Z,G,CP,RKN,WS,EA,...
                         RHOW,HNEU,fstab)

      TAK=TA+TK;
      TSK=TS+TK;
      RKIN=RKINST(RKN,WS,TAK,TSK,Z,G,fstab);
      RHOA=PR/(RA*(TAK));
%   RHOA in kg/m3
      QH=RHOA*(TA-TS)*CP*RKIN;
      ES=SVPI(TS);
      QE=0.622*HNEU/(RA*(TAK))*RKIN*(EA-ES);
      E=-QE/(RHOW*HNEU);
%   E in  m/hr

end