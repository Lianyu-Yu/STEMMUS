%     Routine to correct energy and mass fluxes when
%     numerical overshoots dictate that W was changed in 
%     the calling routine - either because W went negative
%     or due to the liquid fraction being held constant.

function [FM,E,Q,QM,MR,QE,QOTHER]=PREHELP(W1,W,DT,FM1,fac,PS,PRAIN,E,RHOW,HF,Q,QM,...
      QE,HNEU)
%        REAL MR,DT

%  The next statement calculates the value of FM given 
%   the W and w1 values
       FM = (W1-W)/DT*fac-FM1;
%  The next statement calculates the changed MR and E due to the new FM.
%   FM was = PRAIN+PS-MR-E
%   Here the difference is absorbed in MR first and then E if mr < 0.
%   
       MR = max( 0.0 , (PS + PRAIN - FM - E));      
       if MR<=0
           MR=0;
       end
       E = PS + PRAIN - FM - MR;
%   Because melt rate changes the advected energy also changes.  Here
%    advected and melt energy are separated,
       QOTHER = Q + QM - QE;
%    then the part due to melt recalculated
       QM = MR*RHOW*HF;
%    then the part due to evaporation recalculated
       QE = -E*RHOW*HNEU;
%    Energy recombined
       Q = QOTHER - QM + QE;
end