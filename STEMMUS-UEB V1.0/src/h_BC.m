function h_BC
global RHS KL_h Precip Evap NN C4 Trap
global NBCh NBChB BCh BChB hN KT Delt_t DSTOR0 NBChh TIME h_SUR AVAIL0 C4_a
global statevu statevw cump cume cummr outv Ta U HR_a Rns Rn Rnl
global Prain Ps Albedo QH QE E_subl MR QM Q FM tave TSURF QNET QSI QLI vub vw
%%%%%%%%%% Apply the bottom boundary condition called for by NBChB %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if NBChB==1            %-----> Specify matric head at bottom to be ---BChB;
    RHS(1)=BChB;
    C4(1,1)=1;
    RHS(2)=RHS(2)-C4(1,2)*RHS(1);
    C4(1,2)=0;
    C4_a(1)=0;
elseif NBChB==2        %-----> Specify flux at bottom to be ---BChB (Positive upwards);
    RHS(1)=RHS(1)+BChB;
elseif NBChB==3        %-----> NBChB=3,Gravity drainage at bottom--specify flux= hydraulic conductivity;
    RHS(1)=RHS(1)-KL_h(1,1);
end

%%%%%%%%%% Apply the surface boundary condition called for by NBCh  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if NBCh==1             %-----> Specified matric head at surface---equal to hN;
    RHS(NN)=h_SUR(KT);
    C4(NN,1)=1;
    RHS(NN-1)=RHS(NN-1)-C4(NN-1,2)*RHS(NN);
    C4(NN-1,2)=0;
    C4_a(NN-1)=0;
elseif NBCh==2
    if NBChh==1
        RHS(NN)=hN;
        C4(NN,1)=1;
        RHS(NN-1)=RHS(NN-1)-C4(NN-1,2)*RHS(NN);
        C4(NN-1,2)=0;
    else
        RHS(NN)=RHS(NN)-BCh;   %> and a specified matric head (saturation or dryness)was applied;
    end
else
    Evap_Cal;
    if TIME>=576*1800-19*1800 && TIME<=581*1800-19*1800     %TIME>=576*1800 && TIME<=581*1800 15-12-13 0-3:00 p=0.2mm TIME>=466*3600 && TIME<=468*3600
        Precip(KT)=0.02/6/1800;%0.0175/1800;
    elseif TIME>=636*1800-19*1800 && TIME<=641*1800-19*1800   %15-12-14 6-9 p=0.1mm
        Precip(KT)=0.01/6/1800;
    elseif TIME>=828*1800-19*1800 && TIME<=833*1800-19*1800   %15-12-18 6-10 p=9 0.1mm
        Precip(KT)=0.01/6/1800;
    elseif TIME>=858*1800-19*1800 && TIME<=863*1800-19*1800   %15-12-18 21-24 p=0.4mm
        Precip(KT)=0.04/6/1800;
    elseif TIME>=864*1800-19*1800 && TIME<=869*1800-19*1800   %15-12-19 0-3 p=0.1mm
        Precip(KT)=0.01/6/1800;
    elseif TIME>=984*1800-19*1800 && TIME<=989*1800-19*1800   %15-12-21 12-15 p=0.1mm
        Precip(KT)=0.01/6/1800;
    elseif TIME>=1110*1800-19*1800 && TIME<=1115*1800-19*1800   %15-12-24 3-6 p=0.2mm
        Precip(KT)=0.02/6/1800;
    elseif TIME>=1914*1800-19*1800 && TIME<=1919*1800-19*1800   %16-1-9 21-24 p=0.4mm
        Precip(KT)=0.04/6/1800;
    elseif TIME>=1920*1800-19*1800 && TIME<=1925*1800-19*1800   %16-1-10 0-3 p=0.3mm
        Precip(KT)=0.03/6/1800;    % Original 0.3mm
    elseif TIME>=1926*1800-19*1800 && TIME<=1931*1800-19*1800   %16-1-10 3-6 p=0.4mm
        Precip(KT)=0.04/6/1800;
    elseif TIME>=1932*1800-19*1800 && TIME<=1937*1800-19*1800   %16-1-10 6-9 p=0.7mm
        Precip(KT)=0.07/6/1800;
    elseif TIME>=1938*1800-19*1800 && TIME<=1943*1800-19*1800   %16-1-10 9-12 p=0.4mm
        Precip(KT)=0.04/6/1800;
    elseif TIME>=2910*1800-19*1800 && TIME<=2915*1800-19*1800   %16-1-30 15-16 p=0.3mm
        Precip(KT)=0.03/6/1800;
    elseif TIME>=2928*1800-19*1800 && TIME<=2923*1800-19*1800   %16-1-31 0-3 p=0.1mm
        Precip(KT)=0.01/6/1800;
    elseif TIME>=3564*1800-19*1800 && TIME<=3569*1800-19*1800   %16-2-13 6-9 p=0.1mm
        Precip(KT)=0.0/6/1800;
    elseif TIME>=3600*1800-19*1800 && TIME<=3605*1800-19*1800   %16-2-14 0-3 p=0.1mm
        Precip(KT)=0.01/6/1800;
    elseif TIME>=3624*1800-19*1800 && TIME<=3629*1800-19*1800   %16-2-14 12-15 p=1.2mm
        Precip(KT)=0.12/6/1800;
    elseif TIME>=3966*1800-19*1800 && TIME<=3971*1800-19*1800   %16-2-21 15-18 p=0.1mm
        Precip(KT)=0.01/6/1800;%     elseif TIME>=1947*1800 && TIME<=1949*1800   %16-1-10 13.30-16 p=0.3mm
    elseif TIME>=3972*1800-19*1800 && TIME<=3977*1800-19*1800   %16-2-21 18-21 p=0.2mm
        Precip(KT)=0.02/6/1800;
    elseif TIME>=3984*1800-19*1800 && TIME<=3989*1800-19*1800   %16-2-22 0-3 p=0.1mm
        Precip(KT)=0.01/6/1800;
    elseif TIME>=4986*1800-1*1800 && TIME<=4991*1800-1*1800   %16-3-13 21-0 p=0.5mm
        Precip(KT)=0.05/6/1800;%     elseif TIME>=1947*1800 && TIME<=1949*1800   %16-1-10 13.30-16 p=0.3mm
    else
        Precip(KT)=0;
    end
    
    [statevu,statevw,cump,cume,cummr,outv]=snow_Calc(Delt_t,Ta,U,HR_a,Rns,Rnl,KT,Rn,Precip);
    vub(KT)=statevu(KT,1);
    vw(KT)=statevw(KT,1);
    Prain(KT)=outv(KT,1);    % Rainfall(m/h) Part of precipitation modeles as rain
    Ps(KT)=outv(KT,2);       % Snowfall(m/h) Part of precipitation modeles as snow
    Albedo(KT)=outv(KT,3);   % Albedo
    QH(KT)=outv(KT,4);       % Sensible heat flux Qh(KJ/m2/h)
    QE(KT)=outv(KT,5);       % Latent heat flux Qe(KJ/m2/h)
    E_subl(KT)=outv(KT,6);   % Sublimation (m/h)
    MR(KT)=outv(KT,7);       % Melt outflow rate (m/h)
    QM(KT)=outv(KT,8);       % Heat advected by melt outflow(KJ/m2/h)
    Q(KT)=outv(KT,9);        % Total surface energy flux into the snow(KJ/m2/h)
    FM(KT)=outv(KT,10);      % Combined mass fluxes(dW/dt)(m/h)
    tave(KT)=outv(KT,11);    % Snow and underlying soil average temperature T (oC)
    TSURF(KT)=outv(KT,12);   % Snow surface temperature, Ts (oC)
    QNET(KT)=outv(KT,13);    % Net radiation (KJ/m2/h)
    QSI(KT)=outv(KT,14);      % Net shortwave radiation (KJ/m2/h)
    QLI(KT)=outv(KT,15);      % Net longwave radiation (KJ/m2/h)
    if MR(KT)>Ps(KT)
        MR(KT)=Ps(KT);
    end
    if Ps(KT)>0
        AVAIL0=MR(KT)*100/3600+DSTOR0/Delt_t;%MR(KT)*100/3600
    else
        AVAIL0=Precip(KT)+DSTOR0/Delt_t;%MR(KT)*100/3600
    end
    if NBChh==1
        RHS(NN)=hN;
        C4(NN,1)=1;
        RHS(NN-1)=RHS(NN-1)-C4(NN-1,2)*RHS(NN);
        C4(NN-1,2)=0;
        C4_a(NN-1)=0;
    else
        if Ps(KT)>0
            RHS(NN)=RHS(NN)+AVAIL0;
        else
            RHS(NN)=RHS(NN)+AVAIL0-Evap(KT);
        end
    end
end

