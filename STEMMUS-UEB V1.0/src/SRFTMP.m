%%%%%%%%%%%%%%%%%%%%%%%%%%%% SRFTMP () %%%%%%%%%%%%%%%%%%%%%%%%%%*
%     Computes the surface temperature of snow  Tarboton 1996 P16
function [SRFTMP]=SRFTMP(QSI,A,QLI,QPIN,EA,TA,TAVE,TK,PR,RA,CP,RHO,RKN,...
    HNEU,ES,SBC,CS,RS,W,qnetob,IRADFL,WS,Z,G,FC,fstab)
%%This version written on 4/23/97 by Charlie Luce solves directly the
%%energy balance equation using Newtons method - avoiding the linearizations
%%used previously.  The derivative is evaluated numerically over the range
%%ts to fff*ts  where fff = 0.999
fff=0.999;tol=0.0001;
%       TK=273.15;
QSN = QSI*(1.0-A);
%  To ensure all temperatures in kelvin
TAK = TA+TK;
TAVEK  = TAVE+TK;
DENSA  = PR/(RA*TAK);     % Density of Air in kg/m3 if PR in Pa
DENS = RHO;
QP=QPIN;            % store input variable locally without changing global value
if(W<=0 && QP>0),QP=0;end
%   ignore the effect of precip advected
%  energy on the calculation of surface temperature when there is no snow.
%   Without this ridiculously high temperatures can result as the model
%   tries to balance outgoing radiation with precip advected energy.
TS = TAK;                             % first approximation
ER=1.0;
niter = 0;
while(ER>tol && niter<20)
    Tslast = TS;
    [F1] = surfeb(TS,RKN,WS,TAK,Z,G,fstab,QP,DENSA,CP,HNEU,...
        PR,EA,TK,DENS,CS,RS,TAVEK,QSN,QLI,FC,ES,SBC,qnetob,IRADFL);
    [F2] = surfeb(fff*TS,RKN,WS,TAK,Z,G,fstab,QP,DENSA,CP,HNEU,...
        PR,EA,TK,DENS,CS,RS,TAVEK,QSN,QLI,FC,ES,SBC,qnetob,IRADFL);
    TS = TS - ((1-fff) * TS * F1) / (F1 - F2);
    ER = abs(TS - Tslast);
    niter=niter+1;
    if (ER<tol)
        break;
    end
end
TS = TS - TK;
if (W>0 && TS>0)
    SRFTMP = 0;
else
    SRFTMP = TS;
end

end
