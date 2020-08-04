function CondL_h
global Theta_LL Theta_r Theta_s Alpha hh n m Se KL_h Ks MN ML ND DTheta_LLh DTheta_UUh NN RHOL RHOI Theta_U Ratio_ice Imped KT EPCT SAVEKfL_h
global NL J Theta_L h IS KIT TT Thmrlefc CKT CKTN POR KfL_h KfL_T Theta_II Theta_UU T_CRIT L_f g T0 TT_CRIT h_frez hh_frez ISFT Lamda Phi_s SWCC XCAP Theta_cap SFCC Gama_hh hd hm Theta_m KL_h_flm Ks_flm
% PRN The lowest suction head (The maximum value of matric head,considering
% the negative sign before them. The absolute value of which is smallest) at which soil remains saturated.

SFCC=1;
MN=0;
for ML=1:NL
    %     J=IS(ML);
    J=ML;
    for ND=1:2
        MN=ML+ND-1;SAVEKfL_h(ML,ND)=KfL_h(ML,ND);
        if SWCC==1
            if SFCC==1
%          for MN=1:NN
            if abs(hh(MN))>=abs(hd)
%                 Gama_h(MN)=0;
                Gama_hh(MN)=0;
            elseif abs(hh(MN))>=abs(hm)
%                 Gama_h(MN)=1-log(abs(hh(MN)))/log(abs(hm));
%                 Gama_h(MN)=log(abs(hd)/abs(hh(MN)))/log(abs(hd)/abs(hm));
                Gama_hh(MN)=log(abs(hd)/abs(hh(MN)))/log(abs(hd)/abs(hm));
            else
%                 Gama_h(MN)=1;
                Gama_hh(MN)=1;
            end
%          end
% %  Gama_hh(MN)=1;
%          Theta_m(ML)=Gama_hh(MN)*Theta_r(J)+(Theta_s(J)-Gama_hh(MN)*Theta_r(J))*(1+abs(Alpha(J)*(-2))^n(J))^m(J);  %Theta_UU==>Theta_LL   Theta_LL==>Theta_UU
%                 if Theta_m(ML)>=POR(J)
%                     Theta_m(ML)=POR(J);
%                 elseif Theta_m(ML)<=Theta_s(J)
%                     Theta_m(ML)=Theta_s(J);
%                 end
                 Theta_m(ML)=Theta_s(J); %% modified to be consistent with TeC saturated water content, 20190827
                if hh(MN)>=-1e-6
                    Theta_UU(ML,ND)=Theta_s(J);
%                     hh(MN)=-1;
                    DTheta_LLh(ML,ND)=0;%hh_frez(MN)=h_frez(MN);
                    %             Se(ML,ND)=1;
%                     if hh(MN)+hh_frez(MN)<=-1e6
%                         Theta_LL(ML,ND)=Theta_r(J);
%                         DTheta_LLh(ML,ND)=0;
%                         Se(ML,ND)=0;
%                     else
                    if (hh_frez(MN))>=0
                        Theta_LL(ML,ND)=Theta_s(J);
                        DTheta_UUh(ML,ND)=0;%Se(ML,ND)=Theta_LL(ML,ND)/POR(J);
                        Se(ML,ND)=1;
                    else
%                         Theta_LL(ML,ND)=Gama_hh(MN)*Theta_r(J)+(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))/(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^m(J);  %Theta_UU==>Theta_LL   Theta_LL==>Theta_UU
%                         Se(ML,ND)=Theta_LL(ML,ND)/POR(J);
                        if Thmrlefc
                            %                             DTheta_UUh(ML,ND)=(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))*Alpha(J)*n(J)*abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)-1);
%                             if Gama_hh(MN)==1
%                                 Theta_LL(ML,ND)=Gama_hh(MN)*Theta_r(J)+(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))/(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^m(J);  %Theta_UU==>Theta_LL   Theta_LL==>Theta_UU
%                                 DTheta_UUh(ML,ND)=(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))*Alpha(J)*n(J)*abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)-1);
%                                 Se(ML,ND)=Theta_LL(ML,ND)/POR(J);
%                             else
                            if (hh(MN)+hh_frez(MN))<=hd
                                Theta_LL(ML,ND)=0;%Gama_hh(MN)*Theta_r(J)+1e-5;
                                DTheta_UUh(ML,ND)=0;%(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                                Se(ML,ND)=0;
                            else
                                Theta_LL(ML,ND)=Gama_hh(MN)*Theta_r(J)+(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))/(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^m(J);  %Theta_UU==>Theta_LL   Theta_LL==>Theta_UU
                                DTheta_UUh(ML,ND)=(-Theta_r(J))/(abs((hh(MN)+hh_frez(MN)))*log(abs(hd/hm)))*(1-(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)))-Alpha(J)*n(J)*m(J)*(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))*((1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)-1))*(abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^(n(J)-1));  %(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                                %                     DTheta_LLh(ML,ND)=(-Theta_a(J))/(abs(hh(MN))*log(abs(hm)))*(1-(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)))-Alpha(J)*n(J)*m(J)*(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*((1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1))*(abs(Alpha(J)*hh(MN))^(n(J)-1));  %(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                                %                     DTheta_LLh(ML,ND)=(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                                Se(ML,ND)=Theta_LL(ML,ND)/POR(J);
                            end
                            
                        else
                            if abs(hh(MN)-h(MN))<1e-3
                                DTheta_UUh(ML,ND)=(Theta_m(ML)-Theta_r(J))*Alpha(J)*n(J)*abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)-1);
                            else
                                DTheta_UUh(ML,ND)=(Theta_LL(ML,ND)-Theta_L(ML,ND))/(hh(MN)+hh_frez(MN)-h(MN)-h_frez(MN));
                            end
                        end
                        
                    end
%                 elseif hh(MN)<=-1e6
%                     Theta_LL(ML,ND)=Theta_r(J);
%                     Theta_UU(ML,ND)=Theta_r(J);
%                     Theta_II(ML,ND)=0;
% %                     hh(MN)=-1e7;
%                     hh_frez(MN)=0;
%                     DTheta_UUh(ML,ND)=0;
%                     DTheta_LLh(ML,ND)=0;
%                     Se(ML,ND)=0;
                else
%                     Theta_UU(ML,ND)=Gama_hh(MN)*Theta_r(J)+(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))/(1+abs(Alpha(J)*hh(MN))^n(J))^m(J);
                    if Thmrlefc
                        %                         DTheta_LLh(ML,ND)=(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))*Alpha(J)*n(J)*abs(Alpha(J)*(hh(MN)))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*(hh(MN)))^n(J))^(-m(J)-1);
%                         if Gama_hh(MN)==1
%                             Theta_UU(ML,ND)=Gama_hh(MN)*Theta_r(J)+(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))/(1+abs(Alpha(J)*hh(MN))^n(J))^m(J);
%                             DTheta_LLh(ML,ND)=(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
%                         else
                        if Gama_hh(MN)==0
                            Theta_UU(ML,ND)=0;
%                             Theta_LL(ML,ND)=Gama_hh(MN)*Theta_r(J)+1e-5;
                            DTheta_LLh(ML,ND)=0;%(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                            if (hh(MN)+hh_frez(MN))<=hd
                                Theta_LL(ML,ND)=0;
                                DTheta_UUh(ML,ND)=0;%(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                                Se(ML,ND)=0;
                            else
                                Theta_LL(ML,ND)=Gama_hh(MN)*Theta_r(J)+(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))/(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^m(J);  %Theta_UU==>Theta_LL   Theta_LL==>Theta_UU
                                DTheta_UUh(ML,ND)=(-Theta_r(J))/(abs((hh(MN)+hh_frez(MN)))*log(abs(hd/hm)))*(1-(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)))-Alpha(J)*n(J)*m(J)*(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))*((1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)-1))*(abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^(n(J)-1));  %(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                                %                     DTheta_LLh(ML,ND)=(-Theta_a(J))/(abs(hh(MN))*log(abs(hm)))*(1-(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)))-Alpha(J)*n(J)*m(J)*(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*((1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1))*(abs(Alpha(J)*hh(MN))^(n(J)-1));  %(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                                %                     DTheta_LLh(ML,ND)=(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                                Se(ML,ND)=Theta_LL(ML,ND)/POR(J);
                            end
                        else
                            Theta_UU(ML,ND)=Gama_hh(MN)*Theta_r(J)+(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))/(1+abs(Alpha(J)*hh(MN))^n(J))^m(J);
                            DTheta_LLh(ML,ND)=(-Theta_r(J))/(abs(hh(MN))*log(abs(hd/hm)))*(1-(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)))-Alpha(J)*n(J)*m(J)*(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))*((1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1))*(abs(Alpha(J)*hh(MN))^(n(J)-1));  %(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                            %                     DTheta_LLh(ML,ND)=(-Theta_a(J))/(abs(hh(MN))*log(abs(hm)))*(1-(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)))-Alpha(J)*n(J)*m(J)*(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*((1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1))*(abs(Alpha(J)*hh(MN))^(n(J)-1));  %(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                            %                     DTheta_LLh(ML,ND)=(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                            if (hh(MN)+hh_frez(MN))<=hd
                                Theta_LL(ML,ND)=0;
%                                 Theta_LL(ML,ND)=Gama_hh(MN)*Theta_r(J)+1e-5;
                                DTheta_UUh(ML,ND)=0;%(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                                Se(ML,ND)=0;
                            else
                                Theta_LL(ML,ND)=Gama_hh(MN)*Theta_r(J)+(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))/(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^m(J);  %Theta_UU==>Theta_LL   Theta_LL==>Theta_UU
                                DTheta_UUh(ML,ND)=(-Theta_r(J))/(abs((hh(MN)+hh_frez(MN)))*log(abs(hd/hm)))*(1-(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)))-Alpha(J)*n(J)*m(J)*(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))*((1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)-1))*(abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^(n(J)-1));  %(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                                %                     DTheta_LLh(ML,ND)=(-Theta_a(J))/(abs(hh(MN))*log(abs(hm)))*(1-(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)))-Alpha(J)*n(J)*m(J)*(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*((1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1))*(abs(Alpha(J)*hh(MN))^(n(J)-1));  %(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                                %                     DTheta_LLh(ML,ND)=(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                                Se(ML,ND)=Theta_LL(ML,ND)/POR(J);
                            end
                        end
                    else
                        if abs(hh(MN)-h(MN))<1e-3
                            DTheta_LLh(ML,ND)=(Theta_m(ML)-Theta_r(J))*Alpha(J)*n(J)*abs(Alpha(J)*(hh(MN)))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*(hh(MN)))^n(J))^(-m(J)-1);
                        else
                            DTheta_LLh(ML,ND)=(Theta_UU(ML,ND)-Theta_U(ML,ND))/(hh(MN)-h(MN));
                        end
                    end
                    
%                     if hh(MN)+hh_frez(MN)<=-1e6
%                         Theta_LL(ML,ND)=Theta_r(J);
%                         DTheta_LLh(ML,ND)=0;
%                         Se(ML,ND)=0;
%                     else
%                     if (hh(MN)+hh_frez(MN))>=-4
%                         Theta_LL(ML,ND)=Theta_s(J);
%                         DTheta_UUh(ML,ND)=0;%Se(ML,ND)=Theta_LL(ML,ND)/POR(J);
%                         Se(ML,ND)=1;
%                     else
%                         Theta_LL(ML,ND)=Gama_hh(MN)*Theta_r(J)+(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))/(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^m(J);  %Theta_UU==>Theta_LL   Theta_LL==>Theta_UU
%                         if Thmrlefc
% %                             DTheta_UUh(ML,ND)=(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))*Alpha(J)*n(J)*abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)-1);
% %                             if Gama_hh(MN)==1
% %                                 Theta_LL(ML,ND)=Gama_hh(MN)*Theta_r(J)+(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))/(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^m(J);  %Theta_UU==>Theta_LL   Theta_LL==>Theta_UU
% %                                 DTheta_UUh(ML,ND)=(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))*Alpha(J)*n(J)*abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)-1);
% %                                 Se(ML,ND)=Theta_LL(ML,ND)/POR(J);
% %                             else
%                             if (hh(MN)+hh_frez(MN))<=hd
%                                 Theta_LL(ML,ND)=0;
%                                 DTheta_UUh(ML,ND)=0;%(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
%                                 Se(ML,ND)=0;
%                             else
%                                 Theta_LL(ML,ND)=Gama_hh(MN)*Theta_r(J)+(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))/(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^m(J);  %Theta_UU==>Theta_LL   Theta_LL==>Theta_UU
%                                 DTheta_UUh(ML,ND)=(-Theta_r(J))/(abs((hh(MN)+hh_frez(MN)))*log(abs(hd/hm)))*(1-(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)))-Alpha(J)*n(J)*m(J)*(Theta_m(ML)-Gama_hh(MN)*Theta_r(J))*((1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)-1))*(abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^(n(J)-1));  %(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
%                                 %                     DTheta_LLh(ML,ND)=(-Theta_a(J))/(abs(hh(MN))*log(abs(hm)))*(1-(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)))-Alpha(J)*n(J)*m(J)*(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*((1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1))*(abs(Alpha(J)*hh(MN))^(n(J)-1));  %(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
%                                 %                     DTheta_LLh(ML,ND)=(Theta_s(J)-Gama_hh(MN)*Theta_a(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
%                                 Se(ML,ND)=Theta_LL(ML,ND)/POR(J);
%                             end
%                             
%                         else
%                             if abs(hh(MN)-h(MN))<1e-3
%                                 DTheta_UUh(ML,ND)=(Theta_m(ML)-Theta_r(J))*Alpha(J)*n(J)*abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)-1);
%                             else
%                                 %                         DTheta_LLh(ML,ND)=(Theta_LL(ML,ND)-Theta_L(ML,ND))/(hh(MN)-h(MN));
%                                 DTheta_UUh(ML,ND)=(Theta_LL(ML,ND)-Theta_L(ML,ND))/(hh(MN)+hh_frez(MN)-h(MN)-h_frez(MN));
%                             end
%                         end
%                         Se(ML,ND)=Theta_LL(ML,ND)/POR(J);
                        
%                     end
                end
                %%%%%%%%%%%%%%%%%%%% Sin function for ice calculation %%%%%%%%%%%%%%%%%%%%%
            else
                Tf1=273.15+1;Tf2=273.15-3;%XCAP(J)=0.23;
                if hh(MN)>=-1e-6
                    Theta_UU(ML,ND)=Theta_s(J);
                    hh(MN)=-1e-6;
                    DTheta_UUh(ML,ND)=0;
                    if TT(MN)+273.15>Tf1
                        Theta_II(ML,ND)=0;         
                        Theta_LL(ML,ND)=Theta_s(J);

                    elseif TT(MN)+273.15>=Tf2
                        Theta_II(ML,ND)=0.5*(1-sin(pi()*(TT(MN)+273.15-0.5*Tf1-0.5*Tf2)/(Tf1-Tf2)))*XCAP(J);  
                        Theta_LL(ML,ND)=Theta_UU(ML,ND)-Theta_II(ML,ND)*RHOI/RHOL;

                    else
                        Theta_II(ML,ND)=XCAP(J);   
                        Theta_LL(ML,ND)=Theta_UU(ML,ND)-Theta_II(ML,ND)*RHOI/RHOL;

                    end
                    if Theta_LL(ML,ND)<=0.06;
                        Theta_LL(ML,ND)=0.06;
                        DTheta_LLh(ML,ND)=0;
                        Se(ML,ND)=0;
                    else
                        Theta_LL(ML,ND)=Theta_LL(ML,ND);
                        Se(ML,ND)=Theta_LL(ML,ND)/POR(J);
                                DTheta_LLh(ML,ND)=DTheta_UUh(ML,ND);
                    end
                    
                elseif hh(MN)<=-1e7
                    Theta_UU(ML,ND)=Theta_r(J);
                    hh(MN)=-1e7;
                    DTheta_UUh(ML,ND)=0;
                    Theta_II(ML,ND)=0;
                    Theta_LL(ML,ND)=Theta_r(J);
                    Se(ML,ND)=0;
                    DTheta_LLh(ML,ND)=0;
                else
                    Theta_UU(ML,ND)=Theta_r(J)+(Theta_s(J)-Theta_r(J))/(1+abs(Alpha(J)*hh(MN))^n(J))^m(J);
                    
                    if Thmrlefc
                        DTheta_UUh(ML,ND)=(Theta_s(J)-Theta_r(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                    else
                        if abs(hh(MN)-h(MN))<1e-3
                            DTheta_UUh(ML,ND)=(Theta_s(J)-Theta_r(J))*Alpha(J)*n(J)*abs(Alpha(J)*hh(MN))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*hh(MN))^n(J))^(-m(J)-1);
                        else
                            DTheta_UUh(ML,ND)=(Theta_UU(ML,ND)-Theta_U(ML,ND))/(hh(MN)-h(MN));
                        end
                    end
                    if TT(MN)+273.15>Tf1
                        Theta_II(ML,ND)=0;   
                        Theta_LL(ML,ND)=Theta_r(J)+(Theta_s(J)-Theta_r(J))/(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^m(J);  %Theta_UU==>Theta_LL   Theta_LL==>Theta_UU

                    elseif TT(MN)+273.15>=Tf2
                        Theta_II(ML,ND)=0.5*(1-sin(pi()*(TT(MN)+273.15-0.5*Tf1-0.5*Tf2)/(Tf1-Tf2)))*XCAP(J);  
                        Theta_LL(ML,ND)=Theta_UU(ML,ND)-Theta_II(ML,ND)*RHOI/RHOL;%Theta_UU(ML,ND)

                    else
                        Theta_II(ML,ND)=XCAP(J); 
                        Theta_LL(ML,ND)=Theta_UU(ML,ND)-Theta_II(ML,ND)*RHOI/RHOL;%Theta_UU(ML,ND)

                    end
                    if Theta_LL(ML,ND)<=0.06;
                        Theta_LL(ML,ND)=0.06;
                        DTheta_LLh(ML,ND)=0;
                        Se(ML,ND)=0;
                    else
                        Theta_LL(ML,ND)=Theta_LL(ML,ND);
                        Se(ML,ND)=Theta_LL(ML,ND)/POR(J);
                                    DTheta_LLh(ML,ND)=DTheta_UUh(ML,ND);
                    end

                end
            end
        else
            if hh(MN)>=Phi_s(J)
                Theta_UU(ML,ND)=Theta_s(J);
                hh(MN)=Phi_s(J);
                DTheta_UUh(ML,ND)=0;
                if hh(MN)+hh_frez(MN)<=-1e7
                    Theta_LL(ML,ND)=Theta_r(J);
                    DTheta_LLh(ML,ND)=0;
                    Se(ML,ND)=0;
                elseif hh(MN)+hh_frez(MN)>=Phi_s(J)
                    Theta_LL(ML,ND)=Theta_s(J);
                    DTheta_LLh(ML,ND)=0;
                    Se(ML,ND)=1;
                else
                    Theta_LL(ML,ND)=Theta_s(J)*((hh(MN)+hh_frez(MN))/Phi_s(J))^(-1*Lamda(J));%Theta_r(J)+(Theta_s(J)-Theta_r(J))/(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^m(J);  %Theta_UU==>Theta_LL   Theta_LL==>Theta_UU
                    if Thmrlefc
                        DTheta_LLh(ML,ND)=Theta_s(J)/Phi_s(J)*((hh(MN)+hh_frez(MN))/Phi_s(J))^(-1*Lamda(J)-1);%(Theta_s(J)-Theta_r(J))*Alpha(J)*n(J)*abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)-1);
                    else
                        if abs(hh(MN)-h(MN))<1e-3
                            DTheta_LLh(ML,ND)=(Theta_s(J)-Theta_r(J))*Alpha(J)*n(J)*abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)-1);
                        else
                            DTheta_LLh(ML,ND)=(Theta_LL(ML,ND)-Theta_L(ML,ND))/(hh(MN)+hh_frez(MN)-h(MN)-h_frez(MN));
                        end
                    end
                    Se(ML,ND)=Theta_LL(ML,ND)/POR(J);
                end
            elseif hh(MN)<=-1e7
                Theta_LL(ML,ND)=Theta_r(J);
                Theta_UU(ML,ND)=Theta_r(J);
                Theta_II(ML,ND)=0;
                hh(MN)=-1e7;
                hh_frez(MN)=-1e-6;
                DTheta_UUh(ML,ND)=0;
                DTheta_LLh(ML,ND)=0;
                Se(ML,ND)=0;
            else
                Theta_UU(ML,ND)=Theta_s(J)*((hh(MN))/Phi_s(J))^(-1*Lamda(J));%Theta_r(J)+(Theta_s(J)-Theta_r(J))/(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^m(J);  %Theta_UU==>Theta_LL   Theta_LL==>Theta_UU
                if Thmrlefc
                    DTheta_UUh(ML,ND)=Theta_s(J)/Phi_s(J)*((hh(MN))/Phi_s(J))^(-1*Lamda(J)-1);%(Theta_s(J)-Theta_r(J))*Alpha(J)*n(J)*abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)-1);
                else
                    if abs(hh(MN)-h(MN))<1e-3
                        DTheta_UUh(ML,ND)=(Theta_s(J)-Theta_r(J))*Alpha(J)*n(J)*abs(Alpha(J)*(hh(MN)))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*(hh(MN)))^n(J))^(-m(J)-1);
                    else
                        DTheta_UUh(ML,ND)=(Theta_UU(ML,ND)-Theta_U(ML,ND))/(hh(MN)-h(MN));
                    end
                end
                
                if hh(MN)+hh_frez(MN)<=-1e7
                    Theta_LL(ML,ND)=Theta_r(J);
                    DTheta_LLh(ML,ND)=0;
                    Se(ML,ND)=0;
                elseif hh(MN)+hh_frez(MN)>=Phi_s(J)
                    Theta_LL(ML,ND)=Theta_s(J);
                    DTheta_LLh(ML,ND)=0;
                    Se(ML,ND)=1;
                else
                    Theta_LL(ML,ND)=Theta_s(J)*((hh(MN)+hh_frez(MN))/Phi_s(J))^(-1*Lamda(J));%Theta_r(J)+(Theta_s(J)-Theta_r(J))/(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^m(J);  %Theta_UU==>Theta_LL   Theta_LL==>Theta_UU
                    if Thmrlefc
                        DTheta_LLh(ML,ND)=Theta_s(J)/Phi_s(J)*((hh(MN)+hh_frez(MN))/Phi_s(J))^(-1*Lamda(J)-1);%(Theta_s(J)-Theta_r(J))*Alpha(J)*n(J)*abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)-1);
                    else
                        if abs(hh(MN)-h(MN))<1e-3
                            DTheta_LLh(ML,ND)=(Theta_s(J)-Theta_r(J))*Alpha(J)*n(J)*abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^(n(J)-1)*(-m(J))*(1+abs(Alpha(J)*(hh(MN)+hh_frez(MN)))^n(J))^(-m(J)-1);
                        else
                            DTheta_LLh(ML,ND)=(Theta_LL(ML,ND)-Theta_L(ML,ND))/(hh(MN)+hh_frez(MN)-h(MN)-h_frez(MN));
                        end
                    end
                    Se(ML,ND)=Theta_LL(ML,ND)/POR(J);
                    
                end
            end
        end
        %%% 20190911
                    if Theta_UU(ML,ND)<=Gama_hh(MN)*Theta_r(ML);
                        Theta_UU(ML,ND)=Gama_hh(MN)*Theta_r(ML)+1e-5;
                    elseif Theta_UU(ML,ND)>=Theta_s(ML)
                        Theta_UU(ML,ND)=Theta_s(ML);
                    else
                        Theta_UU(ML,ND)=Theta_UU(ML,ND);
                    end
                    if Theta_LL(ML,ND)<=Gama_hh(MN)*Theta_r(ML);
                        Theta_LL(ML,ND)=Gama_hh(MN)*Theta_r(ML)+1e-5;
                    elseif Theta_LL(ML,ND)>=Theta_s(ML)
                        Theta_LL(ML,ND)=Theta_s(ML);
                    else
                        Theta_LL(ML,ND)=Theta_LL(ML,ND);
                    end
        %%%%            
%         if Theta_LL(ML,ND)==Theta_r(J) && Theta_L(ML,ND)>Theta_r(J)
%            Theta_LL(ML,ND)=(Theta_L(ML,ND)+Theta_r(J))/2; 
%            Se(ML,ND)=Theta_L(ML,ND)/POR(J)/2;
% %            Theta_UU(ML,ND)=(Theta_U(ML,ND)+Theta_r(J))/2; 
%         %    Theta_II(ML,ND)=(Theta_I(ML,ND))/2; 
%         %    Delt_t=Delt_t*0.25;
%         %    tS=tS+1;
%         %    return 
%         end
%         if Theta_UU(ML,ND)==Theta_r(J) && Theta_U(ML,ND)>Theta_r(J)
%            Theta_UU(ML,ND)=(Theta_U(ML,ND)+Theta_r(J))/2; 
% %            Se(ML,ND)=Theta_U(ML,ND)/POR(J)/2;
% %            Theta_UU(ML,ND)=(Theta_U(ML,ND)+Theta_r(J))/2; 
%         %    Theta_II(ML,ND)=(Theta_I(ML,ND))/2; 
%         %    Delt_t=Delt_t*0.25;
%         %    tS=tS+1;
%         %    return 
%         end
%         if (hh(MN)+hh_frez(MN))>=-4
%             %Theta_LL(ML,ND)=Theta_s(J);Theta_UU(ML,ND)=Theta_s(J);
%             DTheta_UUh(ML,ND)=0;%Se(ML,ND)=Theta_LL(ML,ND)/POR(J);
%             DTheta_LLh(ML,ND)=0;%Se(ML,ND)=1;
%         end
%         if Gama_hh(MN)==0
%            DTheta_UUh(ML,ND)=0;
%            DTheta_UUh(ML,ND)=0;
%         end
% DTheta_LLh(ML,ND)=DTheta_UUh(ML,ND);
        if Se(ML,ND)>=1
            Se(ML,ND)=1;
        elseif Se(ML,ND)<=0
            Se(ML,ND)=0;
        end
%         if isnan(Se(ML,ND))==1
%             Se(ML,ND)=0;
%         end
        if isnan(Theta_LL(ML,ND))==1
            keyboard%Theta_LL(ML,ND)=0;
        end
        if isnan(Se(ML,ND))==1
            keyboard
        end
        Theta_II(ML,ND)=(Theta_UU(ML,ND)-Theta_LL(ML,ND))*RHOL/RHOI;  % ice water content
        if Theta_UU(ML,ND)~=0
        Ratio_ice(ML,ND)=RHOI*Theta_II(ML,ND)/(RHOL*Theta_UU(ML,ND)); % ice ratio
        else
            Ratio_ice(ML,ND)=0;
        end
%         if KIT
%             MU_W0=2.4152*10^(-4);   %(g.cm^-1.s^-1)
%             MU1=4742.8;                   %(J.mol^-1)
%             MU_WN=MU_W0*exp(MU1/(8.31441*(20+133.3)));
%             MU_W(ML,ND)=MU_W0*exp(MU1/(8.31441*(TT(MN)+133.3)));
%             CKT(MN)=MU_WN/MU_W(ML,ND);
%             if SWCC==1
%                 KL_h(ML,ND)=CKT(MN)*Ks(J)*(Se(ML,ND)^(0.5))*(1-(1-Se(ML,ND)^(1/m(J)))^m(J))^2;
%             else
%                 KL_h(ML,ND)=CKT(MN)*Ks(J)*(Se(ML,ND))^(3+2/Lamda(J));
%             end
%             KfL_h(ML,ND)=KL_h(ML,ND)*10^(-1*Imped(MN)*Ratio_ice(ML,ND));  % hydraulic conductivity for freezing soil
%             KfL_T(ML,ND)=heaviside(TT_CRIT(MN)-(TT(MN)+T0))*L_f*1e4/(g*(T0));   % thermal consider for freezing soil
% %             if KfL_h(ML,ND)<=2e-17
% %                 KfL_h(ML,ND)=2e-17;
% %             end
%         else
%             KL_h(ML,ND)=0;
%             KfL_h(ML,ND)=0;
%             KfL_T(ML,ND)=0;
%         end
        if KIT
%             Sc(ML,ND)=(1+abs(Alpha(J)*he(MN))^n(J))^(-m(J));
            MU_W0=2.4152*10^(-4);   %(g.cm^-1.s^-1)
            MU1=4742.8;                   %(J.mol^-1)
%             
            MU_WN=MU_W0*exp(MU1/(8.31441*(20+133.3)));
%             if TT(MN)<-20
%                 MU_W(ML,ND)=MU_W0*exp(MU1/(8.31441*(TT(MN)+133.3)));
%                 CKT(MN)=MU_WN/MU_W(ML,ND);
%             elseif TT(MN)<=150
%                 MU_W(ML,ND)=MU_W0*exp(MU1/(8.31441*(TT(MN)+133.3)));
%                 CKT(MN)=MU_WN/MU_W(ML,ND);
%                 CKT(MN)=1.5e-6;
%             end
            if TT(MN)<-20
                MU_W(ML,ND)=3.71e-2; %CKT(MN)=0.2688;
            elseif TT(MN)>100
                MU_W(ML,ND)=0.0028;
%                 MU_W(ML,ND)=MU_W0*exp(MU1/(8.31441*(100+133.3)));
                %CKT(MN)=5.5151;   % kg��m^-1��s^-1 --> 10 g.cm^-1.s^-1; J.cm^-2---> kg.m^2.s^-2.cm^-2--> 1e7g.cm^2.s^-2.cm^-2
            else
                MU_W(ML,ND)=MU_W0*exp(MU1/(8.31441*(TT(MN)+133.3)));
                
            end
            
% CKT(MN)=CKTN/(50+2.575*TT(MN));
% if isnan(MU_W(ML,ND))==1
%     MU_W(ML,ND)=MU_WN;
% end
            CKT(MN)=MU_WN/MU_W(ML,ND);
            if isnan(CKT(MN))==1
                CKT(MN)=1;
            end
%             CKT(MN)=1;
            if Se(ML,ND)==0
                KL_h(ML,ND)=0;
            else
                KL_h(ML,ND)=CKT(MN)*Ks(J)*(Se(ML,ND)^(0.5))*(1-(1-Se(ML,ND)^(1/m(J)))^m(J))^2;
            end

            %                 KL_h(ML,ND)=CKT(MN)*Ks(J)*(Se(ML,ND)^(0.5))*((1-(1-(Se(ML,ND)*Sc(ML,ND))^(1/m(J)))^m(J))/(1-(1-(Sc(ML,ND))^(1/m(J)))^m(J)))^2;
            %                 KL_h_flm(ML,ND)=Ks_flm(ML,ND)*(1+RHOL*g*hh(MN)*D/2/sigma);
            CORF=1;
            FILM=1;   % indicator for film flow parameterization; =1, Zhang (2010); =2, Lebeau and Konrad (2010)
            if FILM==1
                % %%%%%%%%% see Zhang (2010)
%                 BP(J)=3.00061397378356e-24; % See page 165 for reference? perfilm.f
                AGR(J)=0.00035;

%                 Coef_Zeta=-1.46789e-5; %see Sutraset -->> fmods_2_2.f
%                 B(J)=BP(J)*(TT(MN)+273.15)^3;
%                 Ks_flm(ML,ND)=CORF*B(J)*(1-POR(J))*(2*AGR(J))^0.5; %m2
%                 if hh(MN)<=0
%                     Kr(ML,ND)=(1+2*AGR(J)*(hh(MN)/100)/Coef_Zeta)^(-1.5);
%                 else
%                     Kr(ML,ND)=1;
%                 end
%                 KL_h_flm(ML,ND)=Ks_flm(ML,ND)*Kr(ML,ND)*1e4; %m2 --> cm2

                %                     % %%%%%%%%% see Zeng (2011) and Zhang (2010)
                %                     BP(J)=3.00061397378356e-24; % See page 165 for reference? perfilm.f
                %                     AGR(J)=0.00035;
                %
%                 Coef_Zeta=-1.46789e-5; %see Sutraset -->> fmods_2_2.f
%                 B(J)=BP(J)*(TT(MN)+273.15)^3;
                RHOW0=1e3;GVA=9.81;
                uw0=2.4152e-5; % Pa s
                u1=4.7428; % kJ mol-1
                R=8.314472; % J mol-1 K-1
                e=78.54;
                e0=8.85e-12;
                B_sig=235.8E-3; 
                Tc=647.15;
                e_sig=1.256;
                c_sig=-0.625;
                kb=1.381e-23; %Boltzmann constant
                Bi=1;Ba=1.602e-19*Bi;
                uw=uw0*exp(u1/R/(TT(MN)+133.3));
                sigma=B_sig*((Tc-(TT(MN)+273.15))/Tc)^e_sig*(1+c_sig*(Tc-(TT(MN)+273.15))/Tc);
                B(J)=2^0.5*pi()^2*RHOW0*GVA/uw*(e*e0/2/sigma)^1.5*(kb*(TT(MN)+273.15)/(Bi*Ba))^3;
                
                
                Coef_f=0.0373*(2*AGR(J))^3.109;
                Ks_flm(ML,ND)=B(J)*(1-POR(J))*(2*AGR(J))^0.5; %m2
                if hh(MN)<-1
                    Kr(ML,ND)=Coef_f*(1+2*AGR(J)*RHOW0*GVA*abs(hh(MN)/100)/2/sigma)^(-1.5);
                else
                    Kr(ML,ND)=1;
                end
                KL_h_flm(ML,ND)=Ks_flm(ML,ND)*Kr(ML,ND)*1e4; %m2 --> cm2
            else
                %%%%%%%%% see sutraset --> perfilm.f based on Lebeau and Konrad (2010)
                %  EFFECTIVE DIAMETER
                ASVL=-6e-20;
                GVA=9.8;
                PERMVAC=8.854e-12; %(VAC)CUM (PERM)EABILITY [PERMFS] [C2J-1M-1 OR F/M OR S4A2/KG/M3]
                ELECTRC=1.602e-19;
                BOTZC=1.381e-23;
                RHOW0=1000;
                PSICM=1000;
                SWM=0.1; %THE SATURATION WHEN SURFACE ROUGHNESS OF THE SOLID GRAIN BECOMES NEGNIGIBLE (-) IT IS EQUIVALENT TO
                %  		SATURATION BEEN 1000M FROM TULLER AND OR (2005)
                RELPW=78.54;
                %                 SWG=1;
                ED=6.0*(1-POR(J))*(-ASVL/(6.0*pi()*RHOW0*GVA*PSICM))^(1.0/3.0)/POR(J)/SWM;
                %  FILM THICKNESS (WARNING) THIS IS NOT THE FULL PART
                DEL=(RELPW*PERMVAC/(2*RHOW0*GVA))^0.50*(pi()*BOTZC*298.15/ELECTRC);
                Ks_flm(ML,ND)=CORF*4.0*(1-POR(J))*DEL^3.0/pi()/ED;

                if hh(MN)<=-1
                    Kr(ML,ND)=(1-Se(ML,ND))*abs(hh(MN)/100)^(-1.50);
                else
                    Kr(ML,ND)=1;
                end
                KL_h_flm(ML,ND)=Ks_flm(ML,ND)*Kr(ML,ND)*1e4; %m2 --> cm2
            end
            if KL_h_flm(ML,ND)<=0
                KL_h_flm(ML,ND)=0;
            elseif KL_h_flm(ML,ND)>=1e-6
                KL_h_flm(ML,ND)=1e-6;
            end
            if KL_h(ML,ND)<=1E-20
                KL_h(ML,ND)=1E-20;
            end
%             if heaviside(TT_CRIT(MN)-(TT(MN)+T0))>0 %&& EPCT(MN)>0%max(EPCT(MN),heaviside(TT_CRIT(MN)-(TT(MN)+T0)))>0
                if Gama_hh(MN)~=1
                    KfL_h(ML,ND)=KL_h(ML,ND)*10^(-1*Imped(MN)*Ratio_ice(ML,ND));%+KL_h_flm(ML,ND);  % hydraulic conductivity for freezing soil
%                                     KL_h(ML,ND)=KL_h(ML,ND)+KL_h_flm(ML,ND);
                else
                    %                 KL_h(ML,ND)=KL_h(ML,ND);
                    KfL_h(ML,ND)=KL_h(ML,ND)*10^(-1*Imped(MN)*Ratio_ice(ML,ND));  % hydraulic conductivity for freezing soil
                end
%             else
% %                 if Gama_hh(MN)~=1
% %                     KL_h(ML,ND)=KL_h(ML,ND)+KL_h_flm(ML,ND);  % hydraulic conductivity for freezing soil
% %                     %                 KL_h(ML,ND)=KL_h(ML,ND)+KL_h_flm(ML,ND);
% %                 else
% %                     %                 KL_h(ML,ND)=KL_h(ML,ND);
% %                     KL_h(ML,ND)=KL_h(ML,ND);  % hydraulic conductivity for freezing soil
% %                 end
%                 KL_h(ML,ND)=KL_h(ML,ND);%+KL_h_flm(ML,ND);
%             end
            if isnan(KL_h(ML,ND))==1
                KL_h(ML,ND)=1e-20;
                keyboard
            end
            if ~isreal(KL_h(ML,ND))
                keyboard
            end
%             if KfL_h(ML,ND)<=1e-20
%                 KfL_h(ML,ND)=1e-20;
%             end
%             KfL_h(ML,ND)=KL_h(ML,ND);
            %             KfL_h(ML,ND)=KL_h(ML,ND)*10^(-1*Imped(MN)*Ratio_ice(ML,ND));  % hydraulic conductivity for freezing soil
            KfL_T(ML,ND)=heaviside(TT_CRIT(MN)-(TT(MN)+T0))*L_f*1e4/(g*(T0));   % thermal consider for freezing soil
%             if KT>120
%                                     if KfL_h(ML,ND)>2e-5
%                                         keyboard
% %                                         KfL_h(ML,ND)=2e-18;
%                                     end
%             end
        else
            KL_h(ML,ND)=0;
            KfL_h(ML,ND)=0;
            KfL_T(ML,ND)=0;

            %                 KL_h(ML,ND)=KL_h(ML,ND);
            %                 1/(Zc-Zd)*(Theta_a(J)+(Theta_s(J)-Theta_a(J))*(1+Alpha(J)^n(J)*exp(n(J)*Zc))^(-m(J)))+n(J)*m(J)*(Theta_s(J)-Theta_a(J))*Alpha(J)^n(J)*exp(n(J)*Zc)/(1+Alpha(J)^n(J)*exp(n(J)*Zc))^(m(J)+1);
        end
    end
end



%%%%%%%%% Unit of KL_h is determined by Ks, which would be given at the%%%%
%%%%%%%%% beginning.Import thing is to keep the unit of matric head hh(MN)
%%%%%%%%% as 'cm'.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%