% Function to calculate Albedo
%    BATS Albedo Model (Dickinson et. al P.21)
function [ALBEDO]=ALBEDO(tausn,coszen,d,aep,abg,AVO,AIRO)

      B = 2.0;
      CS = 0.2;
%     AVO = 0.95
      CN = 0.5;
%     AIRO = 0.65

      FAGE = tausn/(1.0+tausn);

      if (coszen<0.5) 
         FZEN = 1.0/B*((B+1.0)/(1.0+2.0*B*coszen)-1.0);
      else
         FZEN = 0.0;
      end
      AVD = (1.0-CS*FAGE)*AVO;
      AVIS = AVD+0.4*FZEN*(1.0-AVD);
      AIRD = (1.0-CN*FAGE)*AIRO;
      ANIR = AIRD+0.4*FZEN*(1.0-AIRD);
      ALBEDO = (AVIS+ANIR)/2.0;
      if (d<aep)    % need to transition albedo to a bare ground value
        rr=(1.0-d/aep)*exp(-d*0.5/aep);
        ALBEDO=rr*abg+(1.0-rr)*ALBEDO;
      end
end

