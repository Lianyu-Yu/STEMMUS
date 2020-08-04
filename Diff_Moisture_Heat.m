function Diff_Moisture_Heat
global ML MN ND NL nD DTheta_LLh Chh Khh Chg J m n
global Theta_L Theta_LL h hh T TT NN Delt_t RHS SAVE VGm VGn Ts_Min Ts_Max
global C1 C4 C7 DeltZ CTT KTT Lambda_eff c_unsat CHK Theta_r Theta_s Alpha C2 C5 C9
global Srt TT_CRIT L_f RHOI g T0 EPCT hh_frez h_frez Theta_UU Theta_U KfL_h c_L
global KLhBAR KL_T KLTBAR DhDZ DTDZ DTDBAR QL QLT QLH
global CTT_PH DTheta_LLhBAR LHS
Thmrlefct=1;
%%%%%%%%%%%%%%%%%
%   Soil Moisture Part
%%%%%%%%%%%%%%%%%
run  Latent
run  CondT_coeff
for ML=1:NL
    J=ML;
    VGm(J)=m(J);VGn(J)=n(J);
end
Cpcty_Eqn=1;

for MN=1:NN
    if hh(MN)<-1e7
        hh(MN)=-1e7;
    elseif hh(MN)>-1e-6
        hh(MN)=-1e-6;
    end
    
end

% hPARM_Sub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Cpcty_Eqn==1 && Delt_t>1E-5
    MN=0;
    for ML=1:NL
        J=ML;
        for ND=1:2
            MN=ML+ND-1;
            
            if ND==1
                Dth1=(hh(MN)-h(MN))/Delt_t;
                DthU1=(hh(MN)+hh_frez(MN)-h_frez(MN)-h(MN))/Delt_t;
            else
                Dth2=(hh(MN)-h(MN))/Delt_t;
                DthU2=(hh(MN)+hh_frez(MN)-h_frez(MN)-h(MN))/Delt_t;
                
            end
            
            if abs(hh(MN)-h(MN))<1e-6
                DTheta_LLh(ML,ND)=(Theta_s(J)-Theta_r(J))*Alpha(J)*VGn(J)*abs(Alpha(J)*hh(MN))^(VGn(J)-1)*(-VGm(J))*(1+abs(Alpha(J)*hh(MN))^VGn(J))^(-VGm(J)-1);
                DTheta_UUh(ML,ND)=(Theta_s(J)-Theta_r(J))*Alpha(J)*VGn(J)*abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^(VGn(J)-1)*(-VGm(J))*(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^VGn(J))^(-VGm(J)-1);
                
            elseif ND==2
                DtTheta1=(Theta_UU(ML,1)-Theta_U(ML,1))/Delt_t;
                DtTheta2=(Theta_UU(ML,2)-Theta_U(ML,2))/Delt_t;
                DtThetaU1=(Theta_LL(ML,1)-Theta_L(ML,1))/Delt_t;
                DtThetaU2=(Theta_LL(ML,2)-Theta_L(ML,2))/Delt_t;
                if Cpcty_Eqn==1
                    DTheta_LLh(ML,1)=(DtTheta1*(Dth1+5*Dth2)+DtTheta2*(Dth2-Dth1))*(Dth1^2+4*Dth1*Dth2+Dth2^2)^-1;
                    DTheta_LLh(ML,2)=(DtTheta1*(Dth1-Dth2)+DtTheta2*(5*Dth1+Dth2))*(Dth1^2+4*Dth1*Dth2+Dth2^2)^-1;
                    DTheta_UUh(ML,1)=(DtThetaU1*(DthU1+5*DthU2)+DtThetaU2*(DthU2-DthU1))*(DthU1^2+4*DthU1*DthU2+DthU2^2)^-1;
                    DTheta_UUh(ML,2)=(DtThetaU1*(DthU1-DthU2)+DtThetaU2*(5*DthU1+DthU2))*(DthU1^2+4*DthU1*DthU2+DthU2^2)^-1;
                    
                else
                    DTheta_LLh(ML,1)=2*DtTheta1/Dth1-DtTheta2/Dth2;
                    DTheta_LLh(ML,2)=2*DtTheta2/Dth2-DtTheta1/Dth1;
                    DTheta_UUh(ML,1)=2*DtThetaU1/DthU1-DtThetaU2/DthU2;
                    DTheta_UUh(ML,2)=2*DtThetaU2/DthU2-DtThetaU1/DthU1;
                    
                end
            end
        end
    end
else
    MN=0;
    for ML=1:NL
        J=ML;
        for ND=1:2
            MN=ML+ND-1;
            if abs(hh(MN)-h(MN))<1e-3
                DTheta_LLh(ML,ND)=(Theta_s(J)-Theta_r(J))*Alpha(J)*VGn(J)*abs(Alpha(J)*hh(MN))^(VGn(J)-1)*(-VGm(J))*(1+abs(Alpha(J)*hh(MN))^VGn(J))^(-VGm(J)-1);
                DTheta_UUh(ML,ND)=(Theta_s(J)-Theta_r(J))*Alpha(J)*VGn(J)*abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^(VGn(J)-1)*(-VGm(J))*(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^VGn(J))^(-VGm(J)-1);
                
            else
                DTheta_LLh(ML,ND)=(Theta_UU(ML,ND)-Theta_U(ML,ND))/(hh(MN)-h(MN));
                DTheta_UUh(ML,ND)=(Theta_LL(ML,ND)-Theta_L(ML,ND))/(hh(MN)+hh_frez(MN)-h(MN)-+h_frez(MN));
                
            end
        end
    end
    
    if any(isnan(DTheta_LLh))
        keyboard
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for ML=1:NL
    for ND=1:nD
        Chh(ML,ND)=DTheta_LLh(ML,ND);
        Khh(ML,ND)=KfL_h(ML,ND); %
        Chg(ML,ND)=KfL_h(ML,ND);
    end
end

for MN=1:NN              % Clean the space in C1-7 every iteration,otherwise, in *.PARM files,
    for ND=1:2           % C1-7 will be mixed up with pre-storaged data, which will cause extremly crazy for computation, which exactly results in NAN.
        C1(MN,ND)=0;
        C7(MN)=0;
        C4(MN,ND)=0;
        C9(MN)=0; % C9 is the matrix coefficient of root water uptake;
    end
end

for ML=1:NL
    C1(ML,1)=C1(ML,1)+Chh(ML,1)*DeltZ(ML)/2;
    C1(ML+1,1)=C1(ML+1,1)+Chh(ML,2)*DeltZ(ML)/2;%
    
    C4ARG1=(Khh(ML,1)+Khh(ML,2))/(2*DeltZ(ML));%sqrt(Khh(ML,1)*Khh(ML,2))/(DeltZ(ML));%
    C4(ML,1)=C4(ML,1)+C4ARG1;
    C4(ML,2)=C4(ML,2)-C4ARG1;
    C4(ML+1,1)=C4(ML+1,1)+C4ARG1;
    
    C7ARG=(Chg(ML,1)+Chg(ML,2))/2;%sqrt(Chg(ML,1)*Chg(ML,2));%
    C7(ML)=C7(ML)-C7ARG;
    C7(ML+1)=C7(ML+1)+C7ARG;
    
    % Srt, root water uptake;
    C9ARG1=(2*Srt(ML,1)+Srt(ML,2))*DeltZ(ML)/6;%sqrt(Chg(ML,1)*Chg(ML,2));%
    C9ARG2=(Srt(ML,1)+2*Srt(ML,2))*DeltZ(ML)/6;
    C9(ML)=C9(ML)+C9ARG1;
    C9(ML+1)=C9(ML+1)+C9ARG2;
    
end
RHS(1)=-C9(1)-C7(1)+(C1(1,1)*h(1)+C1(1,2)*h(2))/Delt_t;
for ML=2:NL
    RHS(ML)=-C9(ML)-C7(ML)+(C1(ML-1,2)*h(ML-1)+C1(ML,1)*h(ML)+C1(ML,2)*h(ML+1))/Delt_t;
end
RHS(NN)=-C9(NN)-C7(NN)+(C1(NN-1,2)*h(NN-1)+C1(NN,1)*h(NN))/Delt_t;

for MN=1:NN
    for ND=1:2
        C4(MN,ND)=C1(MN,ND)/Delt_t+C4(MN,ND);
    end
end

SAVE(1,1,1)=RHS(1);
SAVE(1,2,1)=C4(1,1);
SAVE(1,3,1)=C4(1,2);
SAVE(2,1,1)=RHS(NN);
SAVE(2,2,1)=C4(NN-1,2);
SAVE(2,3,1)=C4(NN,1);

run h_BC

RHS(1)=RHS(1)/C4(1,1);

for ML=2:NN
    C4(ML,1)=C4(ML,1)-C4(ML-1,2)*C4(ML-1,2)/C4(ML-1,1);
    RHS(ML)=(RHS(ML)-C4(ML-1,2)*RHS(ML-1))/C4(ML,1);
end

for ML=NL:-1:1
    RHS(ML)=RHS(ML)-C4(ML,2)*RHS(ML+1)/C4(ML,1);
end

for MN=1:NN
    CHK(MN)=abs(RHS(MN)-hh(MN));
end

for MN=1:NN
    hh(MN)=RHS(MN);
end
%% 20200324
for ML=1:NL
    DTheta_LLhBAR(ML)=(DTheta_UUh(ML,1)+DTheta_UUh(ML,2))/2;
end
for ML=1:NL
    LHS(ML)= DTheta_LLhBAR(ML)*(hh(ML)-h(ML))/Delt_t;
end
%%
for MN=1:NN
    if hh(MN)<-1e7
        hh(MN)=-1e7;
    elseif hh(MN)>-1e-6
        hh(MN)=-1e-6;
    end
    
end

run h_Bndry_Flux;

%%%%%%%%%%%%%%%%%%
% Heat Transport Part
%%%%%%%%%%%%%%%%%%
if Thmrlefct==1
    for MN=1:NN
        if TT(MN)<=Ts_Min  %0  %
            TT(MN)=Ts_Min; %  0; %
        elseif TT(MN)>=Ts_Max
            TT(MN)=Ts_Max;
        end
    end
    CTg=zeros(NL,nD);
    
    %% calculate QL
    for ML=1:NL
        KLhBAR(ML)=(KfL_h(ML,1)+KfL_h(ML,2))/2;
        KLTBAR(ML)=(KL_T(ML,1)+KL_T(ML,2))/2;
        DhDZ(ML)=(hh(ML+1)-hh(ML))/DeltZ(ML);
        DTDZ(ML)=(TT(ML+1)-TT(ML))/DeltZ(ML);
        DTDBAR(ML)=0;
        DTheta_LLhBAR(ML)=(DTheta_LLh(ML,1)+DTheta_LLh(ML,2))/2;
    end
    
    %%%%%% NOTE: The soil air gas in soil-pore is considered with Xah and XaT terms.(0.0003,volumetric heat capacity)%%%%%%
    MN=0;
    for ML=1:NL
        QL(ML)=-(KLhBAR(ML)*DhDZ(ML)+(KLTBAR(ML)+DTDBAR(ML))*DTDZ(ML)+KLhBAR(ML));
        QLT(ML)=-((KLTBAR(ML)+DTDBAR(ML))*DTDZ(ML));
        QLH(ML)=-(KLhBAR(ML)*DhDZ(ML)+KLhBAR(ML));
    end
    %%
    MN=0;
    for ML=1:NL
        for ND=1:2
            MN=ML+ND-1;
            CTT(ML,ND)=c_unsat(ML,ND);
            KTT(ML,ND)=Lambda_eff(ML,ND);
            if abs(DTheta_LLh(ML,ND)-DTheta_UUh(ML,ND))>0 %max(EPCT(MN),heaviside(TT_CRIT(MN)-(TT(MN)+T0)))>0  %
                CTT_PH(ML,ND)=(10*L_f^2*RHOI/(g*(T0+TT(MN))))*DTheta_UUh(ML,ND)*max(EPCT(MN),heaviside(TT_CRIT(MN)-(TT(MN)+T0)));%*max(EPCT(ML),heaviside(TT_CRIT(MN)-(TT(MN)+T0)))3.85*1e6*DTheta_LLh(ML,ND)*max(EPCT(MN),heaviside(TT_CRIT(MN)-(TT(MN)+T0)));%(10*L_f^2*RHOI/(g*(T0+TT(MN)))-1e4*c_i*(TT(MN))*L_f*RHOI/(g*(T0+TT(MN))))*DTheta_LLh(ML,ND)*max(EPCT(MN),heaviside(TT_CRIT(MN)-(TT(MN)+T0)));
                
                if CTT_PH(ML,ND)<0 %|| Delt_t<1.0e-6
                    CTT_PH(ML,ND)=0;
                else
                    CTT_PH(ML,ND)=CTT_PH(ML,ND);
                end
                CTT(ML,ND)=c_unsat(ML,ND)+CTT_PH(ML,ND);
            else
                CTT(ML,ND)=c_unsat(ML,ND);
            end
            CTg(ML,ND)=-c_L*Srt(ML,ND)*TT(MN);
        end
    end
    
    if any(imag(TT))
        keyboard
    end
    % EnrgyMAT Sub
    for MN=1:NN              % Clean the space in C1-7 every iteration,otherwise, in *.PARM files,
        for ND=1:2           % C1-7 will be mixed up with pre-storaged data, which will cause extremly crazy for computation, which exactly results in NAN.
            C2(MN,ND)=0;
            C5(MN,ND)=0;
        end
    end
    
    for ML=1:NL
        C2(ML,1)=C2(ML,1)+CTT(ML,1)*DeltZ(ML)/2;
        C2(ML+1,1)=C2(ML+1,1)+CTT(ML,2)*DeltZ(ML)/2;
        
        C5ARG1=(KTT(ML,1)+KTT(ML,2))/(2*DeltZ(ML)); %sqrt(KTT(ML,1)*KTT(ML,2))/(DeltZ(ML));%
        C5(ML,1)=C5(ML,1)+C5ARG1;
        C5(ML,2)=C5(ML,2)-C5ARG1;
        C5(ML+1,1)=C5(ML+1,1)+C5ARG1;
        %% RWU root water uptake
        C7ARG=(CTg(ML,1)+CTg(ML,2))/2; %sqrt(CTg(ML,1)*CTg(ML,2));%
        C7(ML)=C7(ML)-C7ARG;
        C7(ML+1)=C7(ML+1)+C7ARG;
    end
    
    % EnrgyEQ_Sub
    RHS(1)=-C7(1)+(C2(1,1)*T(1)+C2(1,2)*T(2))/Delt_t;
    for ML=2:NL
        RHS(ML)=-C7(ML)+(C2(ML-1,2)*T(ML-1)+C2(ML,1)*T(ML)+C2(ML,2)*T(ML+1))/Delt_t;
    end
    RHS(NN)=-C7(NN)+(C2(NN-1,2)*T(NN-1)+C2(NN,1)*T(NN))/Delt_t;
    
    for MN=1:NN
        for ND=1:2
            C5(MN,ND)=C2(MN,ND)/Delt_t+C5(MN,ND);
        end
    end
    
    SAVE(1,1,2)=RHS(1);
    SAVE(1,2,2)=C5(1,1);
    SAVE(1,3,2)=C5(1,2);
    SAVE(2,1,2)=RHS(NN);
    SAVE(2,2,2)=C5(NN-1,2);
    SAVE(2,3,2)=C5(NN,1);
    
    run Enrgy_BC;
    
    RHS(1)=RHS(1)/C5(1,1);
    for ML=2:NN
        C5(ML,1)=C5(ML,1)-C5(ML-1,2)*C5(ML-1,2)/C5(ML-1,1);
        RHS(ML)=(RHS(ML)-C5(ML-1,2)*RHS(ML-1))/C5(ML,1);
    end
    for ML=NL:-1:1
        RHS(ML)=RHS(ML)-C5(ML,2)*RHS(ML+1)/C5(ML,1);
    end
    
    for MN=1:NN
        CHK(MN)=abs(RHS(MN)-TT(MN));
        TT(MN)=RHS(MN);
    end
    
    if any(imag(TT)) || any(isnan(TT))
        keyboard
    end
    
    for MN=1:NN
        if TT(MN)<=Ts_Min  %0 %
            TT(MN)=Ts_Min;  %0; %
        elseif TT(MN)>=Ts_Max
            TT(MN)=Ts_Max;
        end
    end
    run Enrgy_Bndry_Flux;
end
end