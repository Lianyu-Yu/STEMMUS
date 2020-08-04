function Enrgy_BC
global C5 c_L RHOL QMB RHS NN  C5_a L_ts SH Precip L KT
global NBCTB NBCT BCT BCTB DSTOR0 Delt_t T Ts Ta Tbtm
global EVAP Rn Rn_SOIL Resis_a c_a r_a_SOIL TIME LET GSOIL SHH Tbtms SHF Ps tave Albedo MR
%%%%%%%%% Apply the bottom boundary condition called for by NBCTB %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if NBCTB==1
    RHS(1)=Tbtm(KT);
    C5(1,1)=1;
    RHS(2)=RHS(2)-C5(1,2)*RHS(1);
    C5(1,2)=0;
    C5_a(1)=0;
elseif NBCTB==2
    RHS(1)=RHS(1)+BCTB;
else    
    C5(1,1)=C5(1,1)-c_L*RHOL*QMB(KT);
end   

%%%%%%%%%% Apply the surface boundary condition called by NBCT %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if NBCT==1

    RHS(NN)=Ts(KT);%BCT;%30;-MR(KT)
    if (Ps(KT))>0 && RHS(NN)>0
    RHS(NN)=0;
    end
    C5(NN,1)=1;
    RHS(NN-1)=RHS(NN-1)-C5(NN-1,2)*RHS(NN);
    C5(NN-1,2)=0;
    C5_a(NN-1)=0;    
    SHH(KT)=0.1200*c_a*(T(NN)-Ta(KT))/r_a_SOIL(KT);%Resis_a(KT);   % J cm-2 s-1
    SHF(KT)=0.1200*c_a*(T(NN)-Ta(KT))/Resis_a(KT);%Resis_a(KT);   % J cm-2 s-1
elseif NBCT==2
    RHS(NN)=RHS(NN)-BCT;
else
    L_ts(KT)=L(NN); 
    SHH(KT)=0.1200*c_a*(T(NN)-Ta(KT))/r_a_SOIL(KT);%Resis_a(KT);   % J cm-2 s-1
    RHS(NN)=RHS(NN)+100*Rn_SOIL(KT)/3600-RHOL*L_ts(KT)*EVAP(KT)-100*SH(KT)/3600+RHOL*c_L*(Ta(KT)*Precip(KT)+DSTOR0*T(NN)/Delt_t);  
     GSOIL(KT)=100*Rn_SOIL(KT)/3600-RHOL*L_ts(KT)*EVAP(KT)-100*SH(KT)/3600; %Rn_SOIL(KT)-LET(KT)-SH(KT);
end
  
