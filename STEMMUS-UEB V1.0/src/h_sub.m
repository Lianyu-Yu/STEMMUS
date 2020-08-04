function h_sub
global hh MN NN

hPARM;
h_MAT;
h_EQ;
h_BC;
hh_Solve;

for MN=1:NN
    if hh(MN)>=-1e-6
        hh(MN)=-1e-6;
    end
end
if any(isnan(hh)) || any(hh(1:NN)<=-1E12)
    for MN=1:NN
        hh(MN)=hOLD(MN);
    end
end
h_Bndry_Flux;