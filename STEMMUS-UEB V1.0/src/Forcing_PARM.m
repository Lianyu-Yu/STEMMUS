function Forcing_PARM
global Rn TIMEOLD
global Ta Ts U HR_a SH Rns Rnl KT P_Va RHOV_A Rv TopPg h_SUR NBCT
global Ts_msr Tbtm Tb_msr Ta_msr RH_msr Rn_msr WS_msr Pg_msr HourInput Rns_msr SUMTIME NoTime


if TIMEOLD==KT
    Ta(KT)=0;HR_a(KT)=0;Ts(KT)=0;U(KT)=0;SH(KT)=0;Rns(KT)=0;Rnl(KT)=0;Rn(KT)=0;TopPg(KT)=0;h_SUR(KT)=0;
end
if NBCT==1 && KT==1
    Ts(1)=0;
end
HourInput=1;

NoTime(KT)=fix(SUMTIME(KT)/1800);
if NoTime(KT)<=0
    Ta(KT)=Ta_msr(1);
    HR_a(KT)=0.01.*(RH_msr(1));
    U(KT)=100.*(WS_msr(1));
    Rns(KT)=(Rns_msr(1))*8.64/24/100;
    TopPg(KT)=1000.*(Pg_msr(1));
    Ts(KT)=Ts_msr(1);
    Rn(KT)=(Rn_msr(1))*8.64/24/100;
    Tbtm(KT)=Tb_msr(1);
else
    Ta(KT)=Ta_msr(NoTime(KT))+(Ta_msr(NoTime(KT)+1)-Ta_msr(NoTime(KT)))/1800*(SUMTIME(KT)-NoTime(KT)*1800);
    HR_a(KT)=0.01.*(RH_msr(NoTime(KT))+(RH_msr(NoTime(KT)+1)-RH_msr(NoTime(KT)))/1800*(SUMTIME(KT)-NoTime(KT)*1800));
    U(KT)=100.*(WS_msr(NoTime(KT))+(WS_msr(NoTime(KT)+1)-WS_msr(NoTime(KT)))/1800*(SUMTIME(KT)-NoTime(KT)*1800));
    Rns(KT)=(Rns_msr(NoTime(KT))+(Rns_msr(NoTime(KT)+1)-Rns_msr(NoTime(KT)))/1800*(SUMTIME(KT)-NoTime(KT)*1800))*8.64/24/100;
    TopPg(KT)=1000.*(Pg_msr(NoTime(KT))+(Pg_msr(NoTime(KT)+1)-Pg_msr(NoTime(KT)))/1800*(SUMTIME(KT)-NoTime(KT)*1800));
    Ts(KT)=Ts_msr(NoTime(KT))+(Ts_msr(NoTime(KT)+1)-Ts_msr(NoTime(KT)))/1800*(SUMTIME(KT)-NoTime(KT)*1800);
    Rn(KT)=(Rn_msr(NoTime(KT))+(Rn_msr(NoTime(KT)+1)-Rn_msr(NoTime(KT)))/1800*(SUMTIME(KT)-NoTime(KT)*1800))*8.64/24/100;
    Tbtm(KT)=Tb_msr(NoTime(KT))+(Tb_msr(NoTime(KT)+1)-Tb_msr(NoTime(KT)))/1800*(SUMTIME(KT)-NoTime(KT)*1800);
    
end


P_Va(KT)=0.611*exp(17.27*Ta(KT)/(Ta(KT)+237.3))*HR_a(KT);  %The atmospheric vapor pressure (KPa)  (1000Pa=1000.1000.g.100^-1.cm^-1.s^-2)

RHOV_A(KT)=P_Va(KT)*1e4/(Rv*(Ta(KT)+273.15));