% Calculates the vapour pressure at a specified temperature over ice.
% C     using polynomial from Lowe (1977).
function [SVPI]=SVPI(T)
      SVPI=6.109177956 + T * (0.503469897 + T * (0.01886013408 + ...
      T * (0.0004176223716 + T * (5.82472028e-06 + T * ...
        (4.838803174e-08 + T * 1.838826904e-10)))));
      SVPI=SVPI*100;   % convert from mb to Pa
end