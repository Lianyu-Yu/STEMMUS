%%%%%%  TAVG () Calculates the average temperature of snow and interacting soil layer

function [TAVG]=TAVG(UB,W,RHOW,CS,TO,RHOG,DE,CG,HF)

      SNHC = RHOW*W*CS;
      SHC  = RHOG*DE*CG;
      CHC  = SNHC+SHC;

% C      SNHC = Snow heat capacity
% C      SHC  = Soil heat capacity
% C      CHC  = Combined heat capacity  

if UB<=0
    TS=UB/CHC;
else
    AL=UB/(RHOW*HF);
    if W>=AL
        TS=TO;
    else
        TS=(UB-W*RHOW*HF)/CHC;
    end
end
TAVG=TS;
end