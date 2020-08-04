%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%* RKINST() %%%%%%%%%%%%%%%%%%%%%%
      function [RKINST]=RKINST(rkn,ws,ta,ts,z,g,fstab)
%   function to calculate no neutral turbulent transfer coefficient using the 
%   richardson number correction. Tarboton and Luce 1996, UEB, P15
      if ws<=0
        RKINST=0;    %  No wind so no sensible or latent heat fluxes.
      else
        rich=g*(ta-ts)*z/(ws*ws*ta);    % ta must be in K
        if rich>=0 
           RKINST=rkn/(1+10*rich);
        else
          RKINST=rkn*(1-10*rich);
        end
      end
%Linear damping of stability correction through parameter fstab
      RKINST=rkn+fstab*(RKINST-rkn);
      return
      end
