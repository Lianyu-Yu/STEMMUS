%Function to compute gradient of saturated vapour pressure,
%   temperature function over ice
%   Uses Lowe (1977) polynomial

function [DELTA]=DELTA(T)
      a=[0.5030305237,0.0377325502,0.001267995369,...
        2.477563108e-5,3.005693132e-7,2.158542548e-9,...
        7.131097725e-12];
      DELTA=a(1)+T*(a(2)+T*(a(3)+T*(a(4)+T*(a(5)+T*(a(6)+T*a(7))))));
      DELTA=DELTA*100;    % convert from mb to Pa
end