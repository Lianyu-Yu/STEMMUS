%    function to compute surface temperature using Penman/surface
%    resistance analogy and equilibrium approach

function [SNOTMP]= SNOTMP(TSTAR,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,...
RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab)

    QSN = QSI*(1.0-A);    % P9 eq.(10)
    %    To ensure all temperatures in kelvin
    TAK = TA+TK;
    TSTARK = TSTAR+TK;
    TAVEK  = TAVE+TK;
    DENSA  = PR/(RA*TAK);
    %    DENSA in kg/m3 if PR in Pa
    %     CP1 = CP/1000.0
    %    CP1 in kj/kg/oK
    DENS = RHO;
    RKIN=RKINST(RKN,WS,TAK,TSTARK,Z,G,fstab);

    UPPER = QP+(DENSA*CP*TAK)*RKIN...
        -(HNEU*DENSA*0.622*RKIN)/PR*(SVPI(TSTAR)...
        -EA-DELTA(TSTAR)*TSTARK)...
        +DENS*CS*RS*TAVEK;
    DEN = DENS*CS*RS+(DENSA*CP)*RKIN+(DELTA(TSTAR)*...
        HNEU*DENSA*0.622*RKIN)/PR;
    if (IRADFL==0)
        UPPER = UPPER + QSN+QLI+3.0*(1-FC)*ES*SBC*TSTARK^4;
        DEN = DEN+4.0*(1-FC)*ES*SBC*TSTARK^3;
    else
        UPPER = UPPER + qnetob;
    end

    SNOTMP = UPPER/DEN-TK;

end