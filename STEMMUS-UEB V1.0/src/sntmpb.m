% Function to compute surface temperature using bisection

function [sntmpb]=sntmpb(TSTAR,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,...
              RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab)

      QSN = QSI*(1.0-A);
%   To ensure all temperatures in kelvin
      TAK = TA+TK;
      TSTARK = TSTAR+TK;
      TAVEK  = TAVE+TK;
      DENSA  = PR/(RA*TAK);
%   DENSA in kg/m3 if PR in Pa
%    CP1 = CP/1000.0
%   CP1 in kj/kg/oK
      DENS = RHO;
      RKIN=RKINST(RKN,WS,TAK,TSTARK,Z,G,fstab);

      UPPER = QP+(DENSA*CP*TAK)*RKIN...
             -(HNEU*DENSA*0.622*RKIN)/PR*(SVPI(TSTAR)-EA)...
             +DENS*CS*RS*TAVEK;
      DEN = DENS*CS*RS+(DENSA*CP)*RKIN;
      if (IRADFL==0) 
         UPPER = UPPER + QSN+QLI-(1-FC)*ES*SBC*TSTARK^4;
      else
         UPPER = UPPER + qnetob;
      end
      sntmpb=UPPER-DEN*TSTARK;
end