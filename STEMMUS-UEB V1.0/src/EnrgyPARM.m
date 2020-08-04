function EnrgyPARM
global MN ML ND NL hh TT DeltZ P_gg
global CTh CTT CTa KTh KTT KTa VTT VTh VTa CTg Vvh VvT Vaa
global Kcva Kcah KcaT Kcaa Ccah CcaT Ccaa Kaa
global c_a c_L RHOL DRHOVT DRHOVh RHOV Hc RHODA DRHODAz L WW DRHOVhDz DRHOVTDz
global Theta_V Theta_g QL V_A
global KfL_h KL_h KL_T D_Ta Lambda_eff c_unsat D_V Eta D_Vg Xah XaT Xaa DTheta_LLT Soilairefc
global DTheta_LLh DVa_Switch
global Khh KhT Kha KLhBAR KLTBAR DTDBAR DhDZ DTDZ DPgDZ Beta_g DEhBAR DETBAR QV Qa RHOVBAR EtaBAR
global L_f RHOI g c_i KfL_T TT_CRIT T0 EPCT KT DTheta_UUh Srt CTT_PH CTT_LT CTT_g CTT_Lg h_frez hh_frez SFCC XCAP
global QVa QLH QLT DVH DVT QVT QVH DLemdDZ DqLDZ DqVDZ DQaDZ SqVDZ LemdaBAR csBAR CTT_HT csBAR1 csBAR2
%%% KL_h==>KfL_h

for ML=1:NL
    if ~Soilairefc
        KLhBAR(ML)=(KfL_h(ML,1)+KfL_h(ML,2))/2;
        KLTBAR(ML)=(KL_T(ML,1)+KL_T(ML,2))/2;
        DETBAR(ML)=(D_V(ML,1)*Eta(ML,1)+D_V(ML,2)*Eta(ML,2))/2;
        %         DhDZ(ML)=(hh(ML+1)-hh(ML))/DeltZ(ML);
        DhDZ(ML)=(hh(ML+1)+hh_frez(ML+1)-hh(ML)-h_frez(ML))/DeltZ(ML);
        DTDZ(ML)=(TT(ML+1)-TT(ML))/DeltZ(ML);
        DPgDZ(ML)=(P_gg(ML+1)-P_gg(ML))/DeltZ(ML);
    end
    DTDBAR(ML)=(D_Ta(ML,1)+D_Ta(ML,2))/2;
    DEhBAR(ML)=(D_V(ML,1)+D_V(ML,2))/2;
    DRHOVhDz(ML)=(DRHOVh(ML+1)+DRHOVh(ML))/2;
    DRHOVTDz(ML)=(DRHOVT(ML+1)+DRHOVT(ML))/2;
    RHOVBAR(ML)=(RHOV(ML+1)+RHOV(ML))/2;
    EtaBAR(ML)=(Eta(ML,1)+Eta(ML,2))/2;
    LemdaBAR(ML)=(Lambda_eff(ML,1)+Lambda_eff(ML,2))/2; %MODIFY 20200321 LY
    csBAR(ML)=(c_unsat(ML,1)+c_unsat(ML,2))/2; %MODIFY 20200321 LY
end

%%%%%% NOTE: The soil air gas in soil-pore is considered with Xah and XaT terms.(0.0003,volumetric heat capacity)%%%%%%
MN=0;
for ML=1:NL
    for ND=1:2
        MN=ML+ND-1;
        if ~Soilairefc
            QL(ML)=-(KLhBAR(ML)*DhDZ(ML)+(KLTBAR(ML)+DTDBAR(ML))*DTDZ(ML)+KLhBAR(ML));
            QLT(ML)=-((KLTBAR(ML)+DTDBAR(ML))*DTDZ(ML));            
            QLH(ML)=-(KLhBAR(ML)*DhDZ(ML)+KLhBAR(ML));
            Qa(ML)=0;
        else            
            Qa(ML)=-((DEhBAR(ML)+D_Vg(ML))*DRHODAz(ML)-RHODA(ML)*(V_A(ML)+Hc*QL(ML)/RHOL));
        end
        
        if DVa_Switch==1
            QV(ML)=-(DEhBAR(ML)+D_Vg(ML))*DRHOVhDz(ML)*DhDZ(ML)-(DEhBAR(ML)*EtaBAR(ML)+D_Vg(ML))*DRHOVTDz(ML)*DTDZ(ML)+RHOVBAR(ML)*V_A(ML);            
            QVa(ML)=RHOVBAR(ML)*V_A(ML);

        else
            QV(ML)=-(DEhBAR(ML)+D_Vg(ML))*DRHOVhDz(ML)*DhDZ(ML)-(DEhBAR(ML)*EtaBAR(ML)+D_Vg(ML))*DRHOVTDz(ML)*DTDZ(ML);
        end
             DVH(ML)=(DEhBAR(ML))*DRHOVhDz(ML);  
             DVT(ML)=(DEhBAR(ML)*EtaBAR(ML))*DRHOVTDz(ML);
            QVH(ML)=-(DEhBAR(ML)+D_Vg(ML))*DRHOVhDz(ML)*DhDZ(ML);
            QVT(ML)=-(DEhBAR(ML)*EtaBAR(ML)+D_Vg(ML))*DRHOVTDz(ML)*DTDZ(ML);        
        if Soilairefc==1
            Kcah(ML,ND)=c_a*TT(MN)*((D_V(ML,ND)+D_Vg(ML))*Xah(MN)+Hc*RHODA(MN)*KfL_h(ML,ND));
            KcaT(ML,ND)=c_a*TT(MN)*((D_V(ML,ND)+D_Vg(ML))*XaT(MN)+Hc*RHODA(MN)*(KL_T(ML,ND)+D_Ta(ML,ND))); %
            Kcaa(ML,ND)=c_a*TT(MN)*Kaa(ML,ND); %((D_V(ML,ND)+D_Vg(ML))*Xaa(MN)+RHODA(MN)*(Beta_g(ML,ND)+Hc*KL_h(ML,ND)/Gamma_w)); %
            if DVa_Switch==1
                Kcva(ML,ND)=L(MN)*RHOV(MN)*Beta_g(ML,ND);  %(c_V*TT(MN)+L(MN))--->(c_L*TT(MN)+L(MN))
            else
                Kcva(ML,ND)=0;
            end
            Ccah(ML,ND)=c_a*TT(MN)*(-V_A(ML)-Hc*QL(ML)/RHOL)*Xah(MN);
            CcaT(ML,ND)=c_a*TT(MN)*(-V_A(ML)-Hc*QL(ML)/RHOL)*XaT(MN);
            Ccaa(ML,ND)=c_a*TT(MN)*Vaa(ML,ND); %*(-V_A(ML)-Hc*QL(ML)/RHOL)*Xaa(MN); %
        end

        if abs(DTheta_LLh(ML,ND)-DTheta_UUh(ML,ND))~=0%max(EPCT(MN),heaviside(TT_CRIT(MN)-(TT(MN)+T0)))>0
            
            CTT_PH(ML,ND)=(10*L_f^2*RHOI/(g*(T0+TT(MN))))*DTheta_UUh(ML,ND)*max(EPCT(ML),heaviside(TT_CRIT(MN)-(TT(MN)+T0)));%*max(EPCT(ML),heaviside(TT_CRIT(MN)-(TT(MN)+T0)))3.85*1e6*DTheta_LLh(ML,ND)*max(EPCT(MN),heaviside(TT_CRIT(MN)-(TT(MN)+T0)));%(10*L_f^2*RHOI/(g*(T0+TT(MN)))-1e4*c_i*(TT(MN))*L_f*RHOI/(g*(T0+TT(MN))))*DTheta_LLh(ML,ND)*max(EPCT(MN),heaviside(TT_CRIT(MN)-(TT(MN)+T0)));
            CTT_Lg(ML,ND)=(c_L*TT(MN)+L(MN))*Theta_g(ML,ND)*DRHOVT(MN);
            CTT_g(ML,ND)=c_a*TT(MN)*Theta_g(ML,ND)*XaT(MN);
            CTT_LT(ML,ND)=((c_L*TT(MN)-c_i*TT(MN)-WW(ML,ND))*RHOL+((c_L*TT(MN)+L(MN))*RHOV(MN)+c_a*RHODA(MN)*TT(MN))*(RHOL/RHOI-1))*1e4*L_f/(g*(T0+TT(MN)))*DTheta_UUh(ML,ND); %DTheta_LLT(ML,ND)
            CTT_HT(ML,ND)=c_unsat(ML,ND);
            if CTT_PH(ML,ND)<=0
                CTT_PH(ML,ND)=0;
            elseif CTT_PH(ML,ND)>5e3
                CTT_PH(ML,ND)=5e3;
            end
            CTT(ML,ND)=c_unsat(ML,ND)+CTT_Lg(ML,ND)+CTT_g(ML,ND)+CTT_LT(ML,ND)+CTT_PH(ML,ND);
            CTh(ML,ND)=(c_L*TT(MN)+L(MN))*Theta_g(ML,ND)*DRHOVh(MN)+c_a*TT(MN)*Theta_g(ML,ND)*Xah(MN);%;%+c_a*TT(MN)*Theta_g(ML,ND)*Xah(MN)
            CTa(ML,ND)=TT(MN)*Theta_V(ML,ND)*c_a*Xaa(MN);% There is not this term in Milly's work.
        else
            CTT(ML,ND)=c_unsat(ML,ND)+(c_L*TT(MN)+L(MN))*Theta_g(ML,ND)*DRHOVT(MN)+c_a*TT(MN)*Theta_g(ML,ND)*XaT(MN) ...
                +((c_L*TT(MN)-WW(ML,ND))*RHOL-(c_L*TT(MN)+L(MN))*RHOV(MN)-c_a*RHODA(MN)*TT(MN))*DTheta_LLT(ML,ND);
            CTh(ML,ND)=((c_L*TT(MN)-WW(ML,ND))*RHOL-(c_L*TT(MN)+L(MN))*RHOV(MN)-c_a*RHODA(MN)*TT(MN))*DTheta_LLh(ML,ND) ...
                +(c_L*TT(MN)+L(MN))*Theta_g(ML,ND)*DRHOVh(MN)+c_a*TT(MN)*Theta_g(ML,ND)*Xah(MN);%;%+c_a*TT(MN)*Theta_g(ML,ND)*Xah(MN)
            CTa(ML,ND)=TT(MN)*Theta_V(ML,ND)*c_a*Xaa(MN);% There is not this term in Milly's work.
            
            CTT_PH(ML,ND)=0;%(10*L_f^2*RHOI/(g*(T0+TT(MN))))*DTheta_UUh(ML,ND)*max(EPCT(ML),heaviside(TT_CRIT(MN)-(TT(MN)+T0)));%*max(EPCT(ML),heaviside(TT_CRIT(MN)-(TT(MN)+T0)))3.85*1e6*DTheta_LLh(ML,ND)*max(EPCT(MN),heaviside(TT_CRIT(MN)-(TT(MN)+T0)));%(10*L_f^2*RHOI/(g*(T0+TT(MN)))-1e4*c_i*(TT(MN))*L_f*RHOI/(g*(T0+TT(MN))))*DTheta_LLh(ML,ND)*max(EPCT(MN),heaviside(TT_CRIT(MN)-(TT(MN)+T0)));
            CTT_Lg(ML,ND)=(c_L*TT(MN)+L(MN))*Theta_g(ML,ND)*DRHOVT(MN);
            CTT_g(ML,ND)=c_a*TT(MN)*Theta_g(ML,ND)*XaT(MN);
            CTT_LT(ML,ND)=((c_L*TT(MN)-WW(ML,ND))*RHOL-(c_L*TT(MN)+L(MN))*RHOV(MN)-c_a*RHODA(MN)*TT(MN))*DTheta_LLT(ML,ND);%((c_L*TT(MN)-c_i*TT(MN)-WW(ML,ND))*RHOL+((c_L*TT(MN)+L(MN))*RHOV(MN)+c_a*RHODA(MN)*TT(MN))*(RHOL/RHOI-1))*1e4*L_f/(g*(T0+TT(MN)))*DTheta_UUh(ML,ND); %DTheta_LLT(ML,ND)
            CTT_HT(ML,ND)=c_unsat(ML,ND);

        end
        if SFCC==0  %%%%%% ice calculation use Sin function
            if TT(MN)+273.15>Tf1
                CTT_PH(ML,ND)=0;
            elseif TT(MN)+273.15>=Tf2%XCAP(MN)*
                CTT_PH(ML,ND)=L_f*10^(-3)*0.5*cos(pi()*(TT(MN)+273.15-0.5*Tf1-0.5*Tf2)/(Tf1-Tf2))*pi()/(Tf1-Tf2);
            else
                CTT_PH(ML,ND)=0;
            end
            CTT_Lg(ML,ND)=(c_L*TT(MN)+L(MN))*Theta_g(ML,ND)*DRHOVT(MN);
            CTT_g(ML,ND)=c_a*TT(MN)*Theta_g(ML,ND)*XaT(MN);
            CTT_LT(ML,ND)=((c_L*TT(MN)-c_i*TT(MN)-WW(ML,ND))*RHOL+((c_L*TT(MN)+L(MN))*RHOV(MN)+c_a*RHODA(MN)*TT(MN))*(RHOL/RHOI-1))*1e4*L_f/(g*(T0+TT(MN)))*DTheta_UUh(ML,ND); %DTheta_LLT(ML,ND)
            
            CTT(ML,ND)=c_unsat(ML,ND)+CTT_Lg(ML,ND)+CTT_g(ML,ND)+CTT_LT(ML,ND)+CTT_PH(ML,ND);
            CTh(ML,ND)=(c_L*TT(MN)+L(MN))*Theta_g(ML,ND)*DRHOVh(MN)+c_a*TT(MN)*Theta_g(ML,ND)*Xah(MN);%;%+c_a*TT(MN)*Theta_g(ML,ND)*Xah(MN)
            CTa(ML,ND)=TT(MN)*Theta_V(ML,ND)*c_a*Xaa(MN);% There is not this term in Milly's work.
        end
        if CTT(ML,ND)<=0
            CTT(ML,ND)=c_unsat(ML,ND)+(c_L*TT(MN)+L(MN))*Theta_g(ML,ND)*DRHOVT(MN)+c_a*TT(MN)*Theta_g(ML,ND)*XaT(MN) ...
                +((c_L*TT(MN)-WW(ML,ND))*RHOL-(c_L*TT(MN)+L(MN))*RHOV(MN)-c_a*RHODA(MN)*TT(MN))*DTheta_LLT(ML,ND);
        end
        KTh(ML,ND)=L(MN)*(D_V(ML,ND)+D_Vg(ML))*DRHOVh(MN)+c_L*TT(MN)*RHOL*Khh(ML,ND)+Kcah(ML,ND); %; %+Kcah(ML,ND)
        KTT(ML,ND)=Lambda_eff(ML,ND)+c_L*TT(MN)*RHOL*KhT(ML,ND)+KcaT(ML,ND)+L(MN)*(D_V(ML,ND)*Eta(ML,ND)+D_Vg(ML))*DRHOVT(MN);  %;%;  % Revised from: "Lambda_eff(ML,ND)+c_L*TT(MN)*RHOL*KhT(ML,ND);"
        KTa(ML,ND)=Kcva(ML,ND)+Kcaa(ML,ND)+c_L*TT(MN)*RHOL*Kha(ML,ND); % There is not this term in Milly's work.
        
        if DVa_Switch==1
            VTh(ML,ND)=c_L*TT(MN)*RHOL*Vvh(ML,ND)+Ccah(ML,ND)-L(MN)*V_A(ML)*DRHOVh(MN);
            VTT(ML,ND)=c_L*TT(MN)*RHOL*VvT(ML,ND)+CcaT(ML,ND)-L(MN)*V_A(ML)*DRHOVT(MN)-(c_L*(QL(ML)+QV(ML))+c_a*Qa(ML)-2.369*QV(ML));
        else
            VTh(ML,ND)=c_L*TT(MN)*RHOL*Vvh(ML,ND)+Ccah(ML,ND);
            VTT(ML,ND)=c_L*TT(MN)*RHOL*VvT(ML,ND)+CcaT(ML,ND)-(c_L*(QL(ML)+QV(ML))+c_a*Qa(ML)-2.369*QV(ML));
        end
        
        VTa(ML,ND)=Ccaa(ML,ND); %c_a*TT(MN)*Vaa(ML,ND);
        
        CTg(ML,ND)=(c_L*RHOL+c_a*Hc*RHODA(MN))*KfL_h(ML,ND)*TT(MN)-c_L*Srt(ML,ND)*TT(MN); %;;% % Revised from "c_L*T(MN)*KL_h(ML,ND)"
        
    end
end
% for ML=1:NL
%     KTT_LBAR(ML)=(KTT_L(ML,1)+KTT_L(ML,2))/2;
%     KTT_aBAR(ML)=(KTT_a(ML,1)+KTT_a(ML,2))/2;
%     KTT_VBAR(ML)=(KTT_V(ML,1)+KTT_V(ML,2))/2;
% end
for ML=1:NL
    
    DLemdDZ(ML)=(LemdaBAR(ML+1)*DTDZ(ML+1)-LemdaBAR(ML)*DTDZ(ML));%/DeltZ(ML); %20200321, ENRGY_1, conductive heat
    DqLDZ(ML)=-(c_L*TT(ML+1)*QL(ML+1)-c_L*TT(ML)*QL(ML));%/DeltZ(ML);%(KTT_LBAR(ML+1)*DTDZ(ML+1)-KTT_LBAR(ML)*DTDZ(ML))/DeltZ(ML); %          %20200321, ENRGY_2, liquid water convective heat
    DqVDZ(ML)=-(c_L*TT(ML+1)*QV(ML+1)-c_L*TT(ML)*QV(ML));%/DeltZ(ML);%(KTT_VBAR(ML+1)*DTDZ(ML+1)-KTT_VBAR(ML)*DTDZ(ML))/DeltZ(ML); %          %20200321, ENRGY_3, liquid water convective heat
    SqVDZ(ML)=-(L(ML+1)*QV(ML+1)-L(ML)*QV(ML))/DeltZ(ML);                    %20200321, ENRGY_4, VAPOR convective heat
    DQaDZ(ML)=-(c_a*TT(ML+1)*Qa(ML+1)-c_a*TT(ML)*Qa(ML));%/DeltZ(ML);%(KTT_aBAR(ML+1)*DTDZ(ML+1)-KTT_aBAR(ML)*DTDZ(ML))/DeltZ(ML); %          %20200321, ENRGY_3, liquid water convective heat
    csBAR2(ML)=(CTT(ML,1)+CTT(ML,2))/2; %MODIFY 20200321 LY
    csBAR1(ML)=(c_unsat(ML,1)+CTT_PH(ML,1)+c_unsat(ML,2)+CTT_PH(ML,2))/2; %MODIFY 20200321 LY
%     cshBAR(ML)=(CTh(ML,1)+CTh(ML,2))/2; %MODIFY 20200321 LY
%     SHLHSh(ML)=cshBAR(ML)*(hh(ML)-h(ML))/Delt_t;
end
