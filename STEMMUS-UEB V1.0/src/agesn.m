%******************************** AGESN () *****************************
%    Function to calculate Dimensionless age of snow for use in
%    BATS Albedo Model (Dickinson et. al P.21)

function [tausn]=agesn(tausn,DT,PS,TSURF,TK)
TSK=TSURF+TK;  % Express surface temperature in kelvin
R1 = exp(5000.0*(1.0/TK - 1.0/TSK));
R2 = R1^10;
if (R2>1.0)
    R2 = 1.0;
end
R3 = 0.3;
%  Dickinson p 23 gives decay as DT*(R1+R2+R3)/tau0  with
%   tau0 = 10**6 sec.  Here 0.0036 = 3600 s/hr * 10**-6 s**-1
%   since dt is in hours.
tausn = max((tausn+0.0036*(R1+R2+R3)*DT)*...
    (1.0 - 100.0*PS*DT),0);
%       RETURN
end