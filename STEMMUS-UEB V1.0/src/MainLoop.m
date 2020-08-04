clc;
clear all;
close all;
tic;
run Constants
global i tS KT Delt_t TEND TIME MN NN NL ML ND hOLD TOLD h hh T TT P_gOLD P_g P_gg Delt_t0 g
global KIT NIT TimeStep Processing
global SUMTIME hhh TTT P_ggg Theta_LLL DSTOR Thmrlefc CHK Theta_LL Theta_L Theta_UUU Theta_UU Theta_U Theta_III Theta_II Theta_I
global NBCh AVAIL Evap DSTOR0 EXCESS QMT RS BCh hN hSAVE NBChh DSTMAX Soilairefc Trap
global TSAVE IRPT1 IRPT2 AVAIL0 TIMEOLD TIMELAST SRT ALPHA BX alpha_h bx Srt KfL_hh KL_hh Chhh ChTT Khhh KhTT CTTT CTT_PH CTT_LT CTT_g CTT_Lg CCTT_PH CCTT_LT CCTT_g CCTT_Lg C_unsat c_unsat
global QL QL_h QL_T QV Qa KL_h Chh ChT Khh KhT SAVEDSTOR SAVEhh SAVEhhh KfL_h TTT_CRIT TT_CRIT T_CRIT TOLD_CRIT
global h_frez hh_frez hhh_frez hOLD_frez ISFT L_f T0 CTT EPCT EPCTT DTheta_LLh DTheta_LLT DTheta_UUh CKT CKTT DDTheta_LLh DDTheta_LLT DDTheta_UUh Lambda_Eff Lambda_eff EfTCON TTETCON TETCON NBCTB Tbtms Cor_TIME DDhDZ DhDZ DDTDZ DTDZ DDRHOVZ DRHOVZ
global DDEhBAR DEhBAR DDRHOVhDz DRHOVhDz EEtaBAR EtaBAR DD_Vg D_Vg DDRHOVTDz DRHOVTDz KKLhBAR KLhBAR KKLTBAR KLTBAR DDTDBAR DTDBAR QVV QLL CChh
global CChT QVT QVH QVTT QVHH SSa Sa HRA HR QVAA QVa QLH QLT QLHH QLTT DVH DVHH DVT DVTT SSe Se QAA QAa QL_TT QL_HH QLAA QL_a DPgDZ DDPgDZ k_g kk_g VV_A V_A
global Ts_Min Ts_Max
global statevu statevw cump cume cummr outv Ta U HR_a Rns Rn Rnl
global Prain Ps Albedo QH QE E_subl MR QM Q FM tave TSURF QNET QSI QLI vub vw hd
Ts_Min=-50;Ts_Max=80;
CPLD=0; % setting the coupled water and heat transfer; CPLD=1, water and heat transfer are coupled; =0, water and heat transfer are independent
run StartInit;   % Initialize Temperature, Matric potential and soil air pressure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TIMEOLD=0;
TIMELAST=0;
SAVEtS=tS;    
L_f=3.34*1e5; %latent heat of freezing fusion J Kg-1
T0=273.15; % unit K

for i=1:1:1778900                         % Notice here: In this code, the 'i' is strictly used for timestep loop and the arraies index of meteorological forcing data.
    KT=KT+1                          % Counting Number of timesteps
    if KT>1 && Delt_t>(TEND-TIME)
        Delt_t=TEND-TIME;           % If Delt_t is changed due to excessive change of state variables, the judgement of the last time step is excuted.
    end
    TIME=TIME+Delt_t;               % The time elapsed since start of simulation
    TimeStep(KT,1)=Delt_t;
    SUMTIME(KT,1)=TIME;
    Processing=TIME/TEND
    %%%%% Updating the state variables. %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if IRPT1==0 && IRPT2==0
        for MN=1:NN
            
            hOLD_frez(MN)=h_frez(MN);
            h_frez(MN)=hh_frez(MN);
            hhh_frez(MN,KT)=hh_frez(MN);
            TOLD_CRIT(MN)=T_CRIT(MN);
            T_CRIT(MN)=TT_CRIT(MN);
            TTT_CRIT(MN,KT)=TT_CRIT(MN);
            
            hOLD(MN)=h(MN);
            h(MN)=hh(MN);
            hhh(MN,KT)=hh(MN);
            KfL_hh(MN,KT)=KfL_h(MN,2);
            KL_hh(MN,KT)=KL_h(MN,2);
            Chhh(MN,KT)=Chh(MN,2);
            ChTT(MN,KT)=ChT(MN,2);
            Khhh(MN,KT)=Khh(MN,2);
            KhTT(MN,KT)=KhT(MN,2);
            CTTT(MN,KT)=CTT(MN,2);
            EPCTT(MN,KT)=EPCT(MN);
            C_unsat(MN,KT)=c_unsat(MN,2);
            CCTT_PH(MN,KT)=CTT_PH(MN,2);
            CCTT_Lg(MN,KT)=CTT_Lg(MN,2);
            CCTT_g(MN,KT)=CTT_g(MN,2);
            CCTT_LT(MN,KT)=CTT_LT(MN,2);
            Lambda_Eff(MN,KT)=Lambda_eff(MN,2);
            EfTCON(MN,KT)=EfTCON(MN,2);
            TTETCON(MN,KT)=TETCON(MN,2);
            DDhDZ(MN,KT)=DhDZ(MN);
            DDTDZ(MN,KT)=DTDZ(MN);
            DDTDZ(MN,KT)=DTDZ(MN);
            DDRHOVZ(MN,KT)=DRHOVZ(MN);
            if Thmrlefc==1
                TOLD(MN)=T(MN);
                T(MN)=TT(MN);
                TTT(MN,KT)=TT(MN);
            end
            if Soilairefc==1
                P_gOLD(MN)=P_g(MN);
                P_g(MN)=P_gg(MN);
                P_ggg(MN,KT)=P_gg(MN);
            end
            if rwuef==1
                SRT(MN,KT)=Srt(MN,1);
                ALPHA(MN,KT)=alpha_h(MN,1);
                BX(MN,KT)=bx(MN,1);
            end
        end
        DSTOR0=DSTOR;
        
        if KT>1
            run SOIL1
        end
    end
    
    for ML=1:NL
        QVV(ML,KT)=QV(ML);
        QLL(ML,KT)=QL(ML,1);
        DDEhBAR(ML,KT)=DEhBAR(ML);
        DDRHOVhDz(ML,KT)=DRHOVhDz(ML);
        DDhDZ(ML,KT)=DhDZ(ML);
        EEtaBAR(ML,KT)=EtaBAR(ML);
        DD_Vg(ML,KT)=D_Vg(ML);
        DDRHOVTDz(ML,KT)=DRHOVTDz(ML);
        DDTDZ(ML,KT)=DTDZ(ML);
        DDPgDZ(ML,KT)=DPgDZ(ML);
        KKLhBAR(ML,KT)=KLhBAR(ML);
        KKLTBAR(ML,KT)=KLTBAR(ML);
        DDTDBAR(ML,KT)=DTDBAR(ML);
        QVAA(ML,KT)=QVa(ML);
        QAA(ML,KT)=Qa(ML);
        if ~Soilairefc
            QLHH(ML,KT)=QLH(ML);
            QLTT(ML,KT)=QLT(ML);
        else
            QLHH(ML,KT)=QL_h(ML);
            QLTT(ML,KT)=QL_T(ML);
            QLAA(ML,KT)=QL_a(ML);
        end
        DVHH(ML,KT)=DVH(ML);
        DVTT(ML,KT)=DVT(ML);
        SSe(ML,KT)=Se(ML,1);
        SSa(ML,KT)=Sa(ML,1);
        
        QVHH(ML,KT)=QVH(ML);
        QVTT(ML,KT)=QVT(ML);
        DDTheta_LLh(ML,KT)=DTheta_LLh(ML,1);
        DDTheta_LLT(ML,KT)=DTheta_LLT(ML,1);
        CChh(ML,KT)=Chh(ML,1);
        kk_g(ML,KT)=k_g(ML,1);
        VV_A(ML,KT)=V_A(ML);
    end
    if Delt_t~=Delt_t0
        for MN=1:NN
            hh(MN)=h(MN)+(h(MN)-hOLD(MN))*Delt_t/Delt_t0;
            TT(MN)=T(MN)+(T(MN)-TOLD(MN))*Delt_t/Delt_t0;
        end
    end
    hSAVE=hh(NN);
    TSAVE=TT(NN);
    if NBCh==1
        hN=BCh;
        hh(NN)=hN;
        hSAVE=hN;
    elseif NBCh==2
        if NBChh~=2
            if BCh<0
                hN=DSTOR0;
                hh(NN)=hN;
                hSAVE=hN;
            else
                hN=-1e6;
                hh(NN)=hN;
                hSAVE=hN;
            end
        end
    else
        if NBChh~=2
            hN=DSTOR0;
            hh(NN)=hN;
            hSAVE=hN;
        end
    end
    run Forcing_PARM
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for KIT=1:NIT   % Start the iteration procedure in a time step.
        for MN=1:NN
            HH(MN)=hh(MN);
            TT_CRIT(MN)=T0;
            if hh(MN)<=-10^(7)
                TT_CRIT(MN)=-inf;  % unit K
            elseif hh(MN)>=-1e-6
                TT_CRIT(MN)=T0;
            else
                TT_CRIT(MN)=T0+g*T0*hh(MN)/L_f/1e4/2;
            end
            if TT_CRIT(MN)<=T0-6
                TT_CRIT(MN)=T0-6;  % unit K
            elseif TT_CRIT(MN)>=T0
                TT_CRIT(MN)=T0;
            else
                TT_CRIT(MN)=TT_CRIT(MN);
            end
            hh_frez(MN)=L_f*1e4*((TT(MN)+T0)-TT_CRIT(MN))/g/T0*heaviside(TT_CRIT(MN)-(TT(MN)+T0));
            if HH(MN)<=hd
                hh_frez(MN)=0;
            end
            if hh_frez(MN)>=-1e-6
                hh_frez(MN)=0;
            elseif hh_frez(MN)<=-1e6
                hh_frez(MN)=-1e6;
            end
        end
        if IRPT1==1
            if EXCESS<0,NBChh=2;end
        end
        run SOIL2;
        run CondL_T;
        run Density_V;
        run CondL_Tdisp;
        if CPLD==1
            run Latent;
            run Density_DA;
            run CondT_coeff;
            run Condg_k_g;
            run CondV_DE;
            run CondV_DVg;
            
            run h_sub;
        else
            run Diff_Moisture_Heat;
        end
        
        if NBCh==1
            DSTOR=0;
            RS=0;
        elseif NBCh==2
            AVAIL=-BCh;
            EXCESS=(AVAIL+QMT(KT))*Delt_t;
            if abs(EXCESS/Delt_t)<=1e-10,EXCESS=0;end
            DSTOR=min(EXCESS,DSTMAX);
            RS=(EXCESS-DSTOR)/Delt_t;
        else
            AVAIL=AVAIL0-Evap(KT);
            EXCESS=(AVAIL+QMT(KT))*Delt_t;
            if abs(EXCESS/Delt_t)<=1e-10,EXCESS=0;end
            DSTOR=0;
            RS=0;
        end
        
        if CPLD==1
            if Soilairefc==1
                run Air_sub;
            end
            
            if Thmrlefc==1
                run Enrgy_sub;
            end
        end
        
        if max(CHK)<0.001
            break
        end
        hSAVE=hh(NN);
        TSAVE=TT(NN);
    end
    TIMEOLD=KT;
    KIT
    KIT=0;
    run SOIL2;
    if KT>15
        run TimestepCHK
    end
    
    SAVEtS=tS;
    if IRPT1==0 && IRPT2==0
        if KT        % In case last time step is not convergent and needs to be repeated.
            MN=0;
            for ML=1:NL
                for ND=1:2
                    MN=ML+ND-1;
                    Theta_LLL(ML,ND,KT)=Theta_LL(ML,ND);
                    Theta_L(ML,ND)=Theta_LL(ML,ND);
                    Theta_UUU(ML,ND,KT)=Theta_UU(ML,ND);
                    Theta_U(ML,ND)=Theta_UU(ML,ND);
                    Theta_III(ML,ND,KT)=Theta_II(ML,ND);
                    Theta_I(ML,ND)=Theta_II(ML,ND);
                    DDTheta_LLh(ML,KT)=DTheta_LLh(ML,2);
                    DDTheta_UUh(ML,KT)=DTheta_UUh(ML,2);
                end
            end
            run ObservationPoints
        end
        if (TEND-TIME)<1E-3
            for MN=1:NN
                hOLD(MN)=h(MN);
                h(MN)=hh(MN);
                hhh(MN,KT)=hh(MN);
                HRA(MN,KT)=HR(MN);
                if Thmrlefc==1
                    TOLD(MN)=T(MN);
                    T(MN)=TT(MN);
                    TTT(MN,KT)=TT(MN);
                    TOLD_CRIT(MN)=T_CRIT(MN);
                    T_CRIT(MN)=TT_CRIT(MN);
                    TTT_CRIT(MN,KT)=TT_CRIT(MN);
                    hOLD_frez(MN)=h_frez(MN);
                    h_frez(MN)=hh_frez(MN);
                    hhh_frez(MN,KT)=hh_frez(MN);
                end
                if Soilairefc==1
                    P_gOLD(MN)=P_g(MN);
                    P_g(MN)=P_gg(MN);
                    P_ggg(MN,KT)=P_gg(MN);
                end
            end
            break
        end
    end
    if KT>0
        for MN=1:NN
            if CPLD==1
                QL(MN,KT)=QL(MN);
                QL_HH(MN,KT)=QL_h(MN);
                QL_TT(MN,KT)=QL_T(MN);
                QV(MN,KT)=QV(MN);
                SAVEhhh(MN,KT)=SAVEhh(MN);
            end
        end
        SAVEDSTOR(KT)=DSTOR;
    end
end

Computational_Time =toc;
%profile off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('COMPUTATIONAL TIME [s] ')
disp(Computational_Time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CPLD==0
    save('unCPLD21.mat')
elseif CPLD==1 && Soilairefc==1
    save('CPLD_air_snow2.mat')
else
    save('CPLD_noair_snow.mat')
end