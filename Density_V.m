function Density_V
global MN RHOV DRHOVh DRHOVT TT hh HR g Rv RHOV_s DRHOV_sT NN

for MN=1:NN
    %     LHR(MN)=exp((-1e6)*g/(Rv*(TT(MN)+273.15)));
    HR(MN)=exp(hh(MN)*g/(Rv*(TT(MN)+273.15)));
    if HR(MN)<=0.01 || isnan(HR(MN))==1
        HR(MN)=0.01;
    elseif HR(MN)>=1
        HR(MN)=1;
    else
        HR(MN)=exp(hh(MN)*g/(Rv*(TT(MN)+273.15)));
    end
    if TT(MN)<-20
        RHOV_s(MN)=1e-6*exp(31.3716-6014.79/(-20+273.15)-7.92495*0.001*(-20+273.15))/(-20+273.15);
        DRHOV_sT(MN)=RHOV_s(MN)*(6014.79/(-20+273.15)^2-7.92495*0.001)-RHOV_s(MN)/(-20+273.15);
        RHOV(MN)=RHOV_s(MN)*HR(MN);
        DRHOVh(MN)=RHOV_s(MN)*HR(MN)*g/(Rv*(-20+273.15));
        DRHOVT(MN)=RHOV_s(MN)*HR(MN)*(-hh(MN)*g/(Rv*(-20+273.15)^2))+HR(MN)*DRHOV_sT(MN);
    elseif TT(MN)>=100
        RHOV_s(MN)=1e-6*exp(31.3716-6014.79/(100+273.15)-7.92495*0.001*(100+273.15))/(100+273.15);
        DRHOV_sT(MN)=RHOV_s(MN)*(6014.79/(100+273.15)^2-7.92495*0.001)-RHOV_s(MN)/(100+273.15);
        RHOV(MN)=RHOV_s(MN)*HR(MN);
        DRHOVh(MN)=RHOV_s(MN)*HR(MN)*g/(Rv*(100+273.15));
        DRHOVT(MN)=RHOV_s(MN)*HR(MN)*(-hh(MN)*g/(Rv*(100+273.15)^2))+HR(MN)*DRHOV_sT(MN);
    else
        RHOV_s(MN)=1e-6*exp(31.3716-6014.79/(TT(MN)+273.15)-7.92495*0.001*(TT(MN)+273.15))/(TT(MN)+273.15);
        DRHOV_sT(MN)=RHOV_s(MN)*(6014.79/(TT(MN)+273.15)^2-7.92495*0.001)-RHOV_s(MN)/(TT(MN)+273.15);
        RHOV(MN)=RHOV_s(MN)*HR(MN);
        DRHOVh(MN)=RHOV_s(MN)*HR(MN)*g/(Rv*(TT(MN)+273.15));
        DRHOVT(MN)=RHOV_s(MN)*HR(MN)*(-hh(MN)*g/(Rv*(TT(MN)+273.15)^2))+HR(MN)*DRHOV_sT(MN);
    end
    
    if isnan(hh(MN))==1  %hh(MN)<=-1E12 ||
        RHOV(MN)=1e-20;
        DRHOVh(MN)=1e-20;
        DRHOVT(MN)=1e-20;
    end
end