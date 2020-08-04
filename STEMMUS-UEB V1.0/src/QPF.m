% Calculates the heat advected to the snowpack due to rain

function [QPF]=QPF(PR,TA,TO,PS,RHOW,HF,CW,CS)
      if (TA>TO) 
         TRAIN=TA;
         TSNOW=TO;
      else
         TRAIN=TO;
         TSNOW=TA;
      end
      QPF=PR*RHOW*(HF+CW*(TRAIN-TO))+PS*RHOW*CS*(TSNOW-TO);
end