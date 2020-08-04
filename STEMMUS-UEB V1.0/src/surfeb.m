function [surfeb]=surfeb(TS,RKN,WS,TAK,Z,G,fstab,QP,DENSA,CP,HNEU,...
    PR,EA,TK,DENS,CS,RS,TAVEK,QSN,QLI,FC,ES,SBC,qnetob,IRADFL)
%      function to evaluate the surface energy balance for use in solving for
%      surface temperature
%      DGT and C Luce 4/23/97
    [RKIN]=RKINST(RKN,WS,TAK,TS,Z,G,fstab);
    surfeb = QP + RKIN*DENSA*CP*(TAK - TS)...
        +(HNEU*DENSA*0.622*RKIN)/PR*(EA-SVP(TS-TK))...
        -DENS*CS*RS*(TS-TAVEK);
    if IRADFL==0
        surfeb = surfeb + QSN+QLI-(1-FC)*ES*SBC*TS^4;
    else
        surfeb = surfeb + qnetob;
    end

end