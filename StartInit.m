function StartInit

global InitND1 InitND2 InitND3 InitND4 InitND5 BtmT BtmX Btmh% Preset the measured depth to get the initial T, h by interpolation method.
global InitT0 InitT1 InitT2 InitT3 InitT4 InitT5 Dmark
global T MN ML NL NN DeltZ Elmn_Lnth Tot_Depth InitLnth
global InitX0 InitX1 InitX2 InitX3 InitX4 InitX5 Inith0 Inith1 Inith2 Inith3 Inith4 Inith5
global h Theta_s Theta_r m n Alpha Theta_L Theta_LL hh TT P_g P_gg Ks 
global XOLD XWRE NS J POR Thmrlefc IH IS Eqlspace FACc
global porosity SaturatedMC ResidualMC SaturatedK Coefficient_n Coefficient_Alpha
global NBCh NBCT NBCP NBChB NBCTB NBCPB BChB BCTB BCPB BCh BCT BCP BtmPg 
global DSTOR DSTOR0 RS NBChh DSTMAX IRPT1 IRPT2 Soilairefc XK XWILT 
global HCAP TCON SF TCA GA1 GA2 GB1 GB2 S1 S2 HCD TARG1 TARG2 GRAT VPER 
global TERM ZETA0 CON0 PS1 PS2 i KLT_Switch DVT_Switch KaT_Switch
global Kaa_Switch DVa_Switch KLa_Switch
global H1 H2 H3 H4 alpha_h bx LR Lm fr RL0 Srt RL TIME rwuef 


Elmn_Lnth=0;
Dmark=0;

for J=1:NS
    POR(J)=porosity(J);
    Theta_s(J)=SaturatedMC(J);
    Theta_r(J)=ResidualMC(J);
    n(J)=Coefficient_n(J);
    Ks(J)=SaturatedK(J);
    Alpha(J)=Coefficient_Alpha(J);
    m(J)=1-1/n(J);
    XK(J)=0.11; %0.11 This is for silt loam; For sand XK=0.025
    XWILT(J)=Theta_r(J)+(Theta_s(J)-Theta_r(J))/(1+abs(Alpha(J)*(-1.5e4))^n(J))^m(J);
end
    if ~Eqlspace
        Inith0=-(((Theta_s(J)-Theta_r(J))/(InitX0-Theta_r(J)))^(1/m(J))-1)^(1/n(J))/Alpha(J);
        Inith1=-(((Theta_s(J)-Theta_r(J))/(InitX1-Theta_r(J)))^(1/m(J))-1)^(1/n(J))/Alpha(J);
        Inith2=-(((Theta_s(J)-Theta_r(J))/(InitX2-Theta_r(J)))^(1/m(J))-1)^(1/n(J))/Alpha(J);
        Inith3=-(((Theta_s(J)-Theta_r(J))/(InitX3-Theta_r(J)))^(1/m(J))-1)^(1/n(J))/Alpha(J);
        Inith4=-(((Theta_s(J)-Theta_r(J))/(InitX4-Theta_r(J)))^(1/m(J))-1)^(1/n(J))/Alpha(J);
        Inith5=-(((Theta_s(J)-Theta_r(J))/(InitX5-Theta_r(J)))^(1/m(J))-1)^(1/n(J))/Alpha(J);
        Btmh=-(((Theta_s(J)-Theta_r(J))/(BtmX-Theta_r(J)))^(1/m(J))-1)^(1/n(J))/Alpha(J);

        if Btmh==-inf
            Btmh=-1e6;
        end
        
        for ML=1:NL
            Elmn_Lnth=Elmn_Lnth+DeltZ(ML);
            InitLnth(ML)=Tot_Depth-Elmn_Lnth;    
            if abs(InitLnth(ML)-InitND5)<1e-10
                for MN=1:(ML+1)
                    T(MN)=BtmT+(MN-1)*(InitT5-BtmT)/ML;
                    h(MN)=(Btmh+(MN-1)*(Inith5-Btmh)/ML);
                    IS(MN)=1;   %%%%%% Index of soil type %%%%%%%
                    IH(MN)=2;   %%%%%% Index of wetting history of soil which would be assumed as dry at the first with the value of 1 %%%%%%%
                end
                Dmark=ML+2;
            end    
            if abs(InitLnth(ML)-InitND4)<1e-10
                for MN=Dmark:(ML+1)
                    T(MN)=InitT5+(MN-Dmark+1)*(InitT4-InitT5)/(ML+2-Dmark);
                    h(MN)=(Inith5+(MN-Dmark+1)*(Inith4-Inith5)/(ML+2-Dmark));
                    IS(MN-1)=1;
                    IH(MN-1)=2;
                end
                Dmark=ML+2;
            end    
            if abs(InitLnth(ML)-InitND3)<1e-10
                for MN=Dmark:(ML+1)
                    T(MN)=InitT4+(MN-Dmark+1)*(InitT3-InitT4)/(ML+2-Dmark);
                    h(MN)=(Inith4+(MN-Dmark+1)*(Inith3-Inith4)/(ML+2-Dmark));
                    IS(MN-1)=1;
                    IH(MN-1)=2;
                end
                Dmark=ML+2;
            end
            if abs(InitLnth(ML)-InitND2)<1e-10
                for MN=Dmark:(ML+1)
                    T(MN)=InitT3+(MN-Dmark+1)*(InitT2-InitT3)/(ML+2-Dmark); 
                    h(MN)=(Inith3+(MN-Dmark+1)*(Inith2-Inith3)/(ML+2-Dmark));
                    IS(MN-1)=1;
                    IH(MN-1)=2;
                end
                Dmark=ML+2;
            end    
            if abs(InitLnth(ML)-InitND1)<1e-10
                for MN=Dmark:(ML+1)
                    T(MN)=InitT2+(MN-Dmark+1)*(InitT1-InitT2)/(ML+2-Dmark);
                    h(MN)=(Inith2+(MN-Dmark+1)*(Inith1-Inith2)/(ML+2-Dmark));
                    IS(MN-1)=1;
                    IH(MN-1)=2;
                end
                Dmark=ML+2;
            end
            if abs(InitLnth(ML))<1e-10
                for MN=Dmark:(NL+1)
                    T(MN)=InitT1+(MN-Dmark+1)*(InitT0-InitT1)/(NL+2-Dmark);
                    h(MN)=(Inith1+(MN-Dmark+1)*(Inith0-Inith1)/(ML+2-Dmark));
                    IS(MN-1)=1;
                    IH(MN-1)=2;
                end
            end
        end
    else
        for MN=1:NN
            h(MN)=-95;
            T(MN)=22;
            TT(MN)=T(MN);
            IS(MN)=1;
            IH(MN)=2;
        end
    end    
% % H1=-15;H2=-30;H4=-15000;  
%  %       H3=-600;
%     
%     % piecewise linear reduction function 
% MN=0;
% for ML=1:NL
%     for ND=1:2
%         MN=ML+ND-1;
%         if hh(MN) >=H1,
%             alpha_h(ML,ND) = 0;
%         elseif  hh(MN) >=H2,
%             alpha_h(ML,ND) = (H1-hh(MN))/(H1-H2);
%         elseif  hh(MN) >=H3,
%             alpha_h(ML,ND) = 1;
%         elseif  hh(MN) >=H4,
%             alpha_h(ML,ND) = (hh(MN)-H4)/(H3-H4);
%         else
%             alpha_h(ML,ND) = 0;
%         end
%     end
% end
% 
% 
% 
% % root lenth distribution
% Lm=127;
% RL0=1;
% r=(Lm-RL0)/59/86400;
% fr=RL0/(RL0+(Lm-RL0)*exp((-1)*(r*TIME)));
% LR=Lm*fr;
% RL=300;
% Elmn_Lnth=0;
%  for ML=1:NL
%          Elmn_Lnth=Elmn_Lnth+DeltZ(ML);
%          
%             if Elmn_Lnth<RL-LR
%                 bx(ML)=0;
%             elseif Elmn_Lnth>=RL-LR && Elmn_Lnth<RL-0.2*LR
%                 bx(ML)=2.0833*(1-(300-Elmn_Lnth)/LR)/LR;
%             else
%                 bx(ML)=1.66667/LR;
%             end
%         for ND=1:2
%             MN=ML+ND-1;
%             bx(ML,ND)=bx(MN);
%         end
%  end
 
for MN=1:NN
    hh(MN)=h(MN);
    if Thmrlefc==1         
        TT(MN)=T(MN);
    end
    if Soilairefc==1
        P_g(MN)=94197.850*10;
        P_gg(MN)=P_g(MN); 
    end
    if MN<NN
    XWRE(MN,1)=0;
    XWRE(MN,2)=0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HCAP(1)=0.998*4.182;HCAP(2)=0.0003*4.182;HCAP(3)=2.66;HCAP(4)=2.66;HCAP(5)=1.3;% ZENG origial HCAP(3)=0.46*4.182;HCAP(4)=0.46*4.182;HCAP(5)=0.6*4.182;    % J cm^-3 Cels^-1  /  g.cm-3---> J g-1 Cels-1;                     %
TCON(1)=1.37e-3*4.182;TCON(2)=6e-5*4.182;TCON(3)=8.8e-2;TCON(4)=2.9e-2;TCON(5)=2.5e-3;% ZENG origial TCON(3)=2.1e-2*4.182;TCON(4)=7e-3*4.182;TCON(5)=6e-4*4.182; % J cm^-1 s^-1 Cels^-1;                %
SF(1)=0;SF(2)=0;SF(3)=0.125;SF(4)=0.125;SF(5)=0.5;                                                                                          %
TCA=6e-5*4.182;GA1=0.035;GA2=0.013;                                                                                                           %
VPER(1)=0.41;VPER(2)=0.05;VPER(3)=0.05;% for sand VPER(1)=0.65;VPER(2)=0;VPER(3)=0;   %  For Silt Loam; % VPER(1)=0.16;VPER(2)=0.33;VPER(3)=0.05;  %
                                                                                                                                                %
%%%%% Perform initial thermal calculations for each soil type. %%%%                                                                             %        
for J=1:NS   %--------------> Sum over all phases of dry porous media to find the dry heat capacity                                             %
    S1=POR(J)*TCA;  %-------> and the sums in the dry thermal conductivity;                                                                     %
    S2=POR(J);                                                                                                                                  %
    HCD(J)=0;                                                                                                                                   %
    for i=3:5                                                                                                                                   %
        TARG1=TCON(i)/TCA-1;                                                                                                                    %
        GRAT=0.667/(1+TARG1*SF(i))+0.333/(1+TARG1*(1-2*SF(i)));                                                                                 %
        S1=S1+GRAT*TCON(i)*VPER(i-2);                                                                                                           %
        S2=S2+GRAT*VPER(i-2);                                                                                                                   %
        HCD(J)=HCD(J)+HCAP(i)*VPER(i-2);                                                                                                        %
    end                                                                                                                                         %
    ZETA0(J)=1/S2;                                                                                                                              %
    CON0(J)=1.25*S1/S2;                                                                                                                         %
    PS1(J)=0;                                                                                                                                   %
    PS2(J)=0;                                                                                                                                   %
    for i=3:5                                                                                                                                   %
        TARG2=TCON(i)/TCON(1)-1;                                                                                                                %
        GRAT=0.667/(1+TARG2*SF(i))+0.333/(1+TARG2*(1-2*SF(i)));                                                                                 %
        TERM=GRAT*VPER(i-2);                                                                                                                    %
        PS1(J)=PS1(J)+TERM*TCON(i);                                                                                                             %
        PS2(J)=PS2(J)+TERM;                                                                                                                     %
    end                                                                                                                                         %
    GB1(J)=0.298/POR(J);                                                                                                                        %
    GB2(J)=(GA1-GA2)/XWILT(J)+GB1(J);                                                                                                           %
end                                                                                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% According to hh value get the Theta_LL
run SOIL2;   % For calculating Theta_LL,used in first Balance calculation.

for ML=1:NL
        Theta_L(ML,1)=Theta_LL(ML,1);
        Theta_L(ML,2)=Theta_LL(ML,2);
        XOLD(ML)=(Theta_L(ML,1)+Theta_L(ML,2))/2;
end
% Using the initial condition to get the initial balance
% information---Initial heat storage and initial moisture storage.
KLT_Switch=1;
DVT_Switch=1;
if Soilairefc
    KaT_Switch=1;
    Kaa_Switch=1;
    DVa_Switch=1;
    KLa_Switch=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% The boundary condition information settings.%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IRPT1=0;
IRPT2=0;
NBCh=3;      % Moisture Surface B.C.: 1 --> Specified matric head(BCh); 2 --> Specified flux(BCh); 3 --> Atmospheric forcing;
BCh=-20/3600;
NBChB=2;    % Moisture Bottom B.C.: 1 --> Specified matric head (BChB); 2 --> Specified flux(BChB); 3 --> Zero matric head gradient (Gravitiy drainage);
BChB=0; 
if Thmrlefc==1
    NBCT=1;  % Energy Surface B.C.: 1 --> Specified temperature (BCT); 2 --> Specified heat flux (BCT); 3 --> Atmospheric forcing;
    BCT=29.75510204;  % surface temperature
    NBCTB=1;% Energy Bottom B.C.: 1 --> Specified temperature (BCTB); 2 --> Specified heat flux (BCTB); 3 --> Zero temperature gradient;
    BCTB=19.8; 
end
if Soilairefc==1
    NBCP=3; % Soil air pressure B.C.: 1 --> Ponded infiltration caused a specified pressure value; 
                % 2 --> The soil air pressure is allowed to escape after beyond the threshold value;
                % 3 --> The atmospheric forcing;
    BCP=0;  
    NBCPB=2;  % Soil air Bottom B.C.: 1 --> Bounded bottom with specified air pressure; 2 --> Soil air is allowed to escape from bottom;
    BCPB=0;  
end

if NBCh~=1
    NBChh=2;                    % Assume the NBChh=2 firstly;
end

FACc=0;                         % Used in MeteoDataCHG for check is FAC changed?
BtmPg=94197.850*10;     % Atmospheric pressure at the bottom (Pa), set fixed
                                     % with the value of mean atmospheric pressure;
DSTOR=0;                        % Depth of depression storage at end of current time step;
DSTOR0=DSTOR;              % Dept of depression storage at start of current time step;
RS=0;                             % Rate of surface runoff;
DSTMAX=0;                     % Depression storage capacity;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



