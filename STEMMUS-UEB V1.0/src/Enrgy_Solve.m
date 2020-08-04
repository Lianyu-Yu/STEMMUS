function Enrgy_Solve
global C5 TT RHS ML NN NL CHK C5_a csBAR Delt_t T HC DeltZ HC1 csBAR1 HC2 csBAR2

RHS(1)=RHS(1)/C5(1,1);

for ML=2:NN
    C5(ML,1)=C5(ML,1)-C5_a(ML-1)*C5(ML-1,2)/C5(ML-1,1);
    RHS(ML)=(RHS(ML)-C5_a(ML-1)*RHS(ML-1))/C5(ML,1);
end

for ML=NL:-1:1
    RHS(ML)=RHS(ML)-C5(ML,2)*RHS(ML+1)/C5(ML,1);
end
for ML=1:NL
    HC(ML)=csBAR(ML)*(TT(ML)-T(ML))/Delt_t*DeltZ(ML);  %20200321 HEAT CONTENT, LHS
    HC1(ML)=csBAR1(ML)*(TT(ML)-T(ML))/Delt_t*DeltZ(ML);  %20200321 HEAT CONTENT, LHS
    HC2(ML)=csBAR2(ML)*(TT(ML)-T(ML))/Delt_t*DeltZ(ML);  %20200321 HEAT CONTENT, LHS

end
for MN=1:NN
    CHK(MN)=abs(RHS(MN)-TT(MN)); %abs((RHS(MN)-TT(MN))/TT(MN)); %
    TT(MN)=RHS(MN);
end
