% Calculates the vapour pressure at a specified temperature over water
% C     using polynomial from Lowe (1977).
function [SVPW]= SVPW(T)
      SVPW=6.107799961 + T * (0.4436518521 + T * (0.01428945805 +... 
      T * (0.0002650648471 + T * (3.031240936e-06 + T * ...
      (2.034080948e-08 + T * 6.136820929e-11)))));
      SVPW=SVPW*100;   % convert from mb to Pa
end