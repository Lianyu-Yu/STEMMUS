%C**************************** atf () ****************************
%    to get the atmospheric transmissivity using the Bristow and Campbell
%    (1984) approach

function [atff]=atf(trange,month,dtbar,a,c)
      b=0.036*exp(-0.154*dtbar(month));
      atff=a*(1-exp(-b*trange^c));
end