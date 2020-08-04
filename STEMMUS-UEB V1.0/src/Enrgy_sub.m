function Enrgy_sub
global TT MN NN TOLD


EnrgyPARM;
Enrgy_MAT;
Enrgy_EQ;
Enrgy_BC;
Enrgy_Solve;
if any(isnan(TT)) || any(abs(TT(1:NN))>=100)
    for MN=1:NN
        TT(MN)=TOLD(MN);
    end
else
    for MN=1:NN
        TT(MN)=TT(MN);
    end
end

Enrgy_Bndry_Flux;