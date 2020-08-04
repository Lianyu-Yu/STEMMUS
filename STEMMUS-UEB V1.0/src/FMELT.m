% Calculates the melt rate and melt outflow

function [FMELT]=FMELT(UB,RHOW,W,HF,LC,RID,KS,PRAIN)
    UU=0.0;
    if (UB<0)
        FMELT=0.0;
    elseif (W<=0.0)
        FMELT=PRAIN;
        if (PRAIN<=0)
            FMELT=0.0;
        end
    else
        UU=UB/(RHOW*W*HF);
        %                           liquid fraction
        if(UU>0.99),UU=0.99;end%end
        %                            TO ENSURE HIGH MELT RATE WHEN ALL LIQUID

        if((UU/(1-UU))<=LC)
            SS=0.0;
        else
            SS=(UU/((1- UU)*RID)-LC/RID)/(1-LC/RID);
        end
    FMELT=KS*SS^3;
    end
    if (FMELT<0.0)
        FMELT=0.0;
    end
end