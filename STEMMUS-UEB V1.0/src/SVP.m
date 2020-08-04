% Calculates the vapour pressure at a specified temperature over water or ice
% c     depending upon temperature.  Temperature is celsius here.
function [SVP]=SVP(T)
    if(T >= 0)
        SVP=SVPW(T);
    else 
        SVP=SVPI(T);
    end
end