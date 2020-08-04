function Dtrmn_Z
%  The determination of the element length
global Elmn_Lnth ML DeltZ NL Tot_Depth DeltZ_R MML 

Elmn_Lnth=0;

for ML=1:3
    DeltZ_R(ML)=0.1;
end
    DeltZ_R(4)=0.2;
    DeltZ_R(5)=0.5;
for ML=6:11
    DeltZ_R(ML)=1;
end
for ML=12:13
    DeltZ_R(ML)=1.5;
end
for ML=14:23
    DeltZ_R(ML)=2;
end
for ML=24:27
    DeltZ_R(ML)=2.5;
end
for ML=28:31
    DeltZ_R(ML)=5;
end
% Sum of element lengths and compared to the total lenght, so that judge 
% can be made to determine the length of rest elements.

for ML=1:31
    Elmn_Lnth=Elmn_Lnth+DeltZ_R(ML);
end

% If the total sum of element lenth is over the predefined depth, stop the
% for loop, make the ML, at which the element lenth sumtion is over defined
% depth, to be new NL.
for ML=32:33
    DeltZ_R(ML)=10;
    Elmn_Lnth=Elmn_Lnth+DeltZ_R(ML);
    if Elmn_Lnth>Tot_Depth
        DeltZ_R(ML)=Tot_Depth-Elmn_Lnth+DeltZ_R(ML);
        NL=ML;

        for ML=1:NL
            MML=NL-ML+1;
            DeltZ(ML)=DeltZ_R(MML);
        end        
        return
    elseif Elmn_Lnth<Tot_Depth && ML==NL
        NL=NL+NL*2;
    end
end

for ML=34:34
    DeltZ_R(ML)=20;
    Elmn_Lnth=Elmn_Lnth+DeltZ_R(ML);
    if Elmn_Lnth>Tot_Depth
        DeltZ_R(ML)=Tot_Depth-Elmn_Lnth+DeltZ_R(ML);
        NL=ML;

        for ML=1:NL
            MML=NL-ML+1;
            DeltZ(ML)=DeltZ_R(MML);
        end        
        return
    elseif Elmn_Lnth<Tot_Depth && ML==NL
        NL=NL+NL*2;
    end
end
    
for ML=35:NL
    DeltZ_R(ML)=20;
    Elmn_Lnth=Elmn_Lnth+DeltZ_R(ML);
    if Elmn_Lnth>=Tot_Depth
        DeltZ_R(ML)=Tot_Depth-Elmn_Lnth+DeltZ_R(ML);
        NL=ML;
        
        for ML=1:NL
            MML=NL-ML+1;
            DeltZ(ML)=DeltZ_R(MML);
        end
        return
    end
end





        