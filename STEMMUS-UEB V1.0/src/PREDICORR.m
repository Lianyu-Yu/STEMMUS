%    Predictor-corrector scheme to update the state variables,
%    U and W for each time step

function [QH,QE,E,MR,QM,Q,FM,TSURF,QNET,W,UB]=PREDICORR(DT,UB,W,A,TA,PRAIN,PS,WS,RH,QSI,...
    QLI,IRADFL,RKN,QNETOB,CDH,CDE,RID,param,sitev)

Wtol=0.025; Utol=2000;
ncall= 0;
ncall=ncall+1;
fac=0.5;
% Parameters
HF =param(7);    % Heat of fusion (333.5 KJ/kg)
RHOW=param(19);  % Density of Water (1000 kg/m^3)
HNEU=param(8);   % Heat of vaporization of ice. (2834 kJ/kg)
[FM,Q,QM,MR,QE,E,TSURF,QH,QNET]=QFM(UB,W,A,TA,PRAIN,PS,WS,RH,QSI,QLI,RKN,IRADFL,...
    QNETOB,CDH,CDE,RID,param,sitev);
W1 = W + DT*FM;
if (W1<0.0)
    W1=0.0;
    [FM,E,Q,QM,MR,QE,QOTHER]=PREHELP(W1,W,DT,0,1,PS,PRAIN,E,RHOW,HF,Q,QM,...
        QE,HNEU);
end
UB1 = UB + DT*Q;
Q1 = Q;
FM1 = FM;
%  save values so that they can be averaged for output
QH1=QH;
QE1=QE;
E1=E;
MR1=MR;
QM1=QM;
TSURF1=TSURF;
QNET1=QNET;

[FM,Q,QM,MR,QE,E,TSURF,QH,QNET]=QFM(UB1,W1,A,TA,PRAIN,PS,WS,RH,QSI,QLI,RKN,IRADFL,...
    QNETOB,CDH,CDE,RID,param,sitev);
W2 = W + DT/2.0*(FM1 + FM);
if (W2<0.0)
    W2=0.0;
    [FM,E,Q,QM,MR,QE,QOTHER]=PREHELP(W2,W,DT,FM1,2,PS,PRAIN,E,RHOW,HF,Q,QM,...
        QE,HNEU);
    
end
UB2 = UB + DT/2.0*(Q1 + Q);
%  iterate to convergence to enhance stability
niter=1;
imax=5;
while((abs(W2-W1)>Wtol || abs(UB2-UB1)>Utol) && ...
        (niter < imax))
    W1=W2;
    UB1=UB2;
    [FM,Q,QM,MR,QE,E,TSURF,QH,QNET]=QFM(UB1,W1,A,TA,PRAIN,PS,WS,RH,QSI,QLI,RKN,IRADFL,...
        QNETOB,CDH,CDE,RID,param,sitev);
    %   corrector again
    W2 = W + DT/2.0*(FM1 + FM);
    if (W2<0.0)
        W2=0.0;
        [FM,E,Q,QM,MR,QE,QOTHER]=PREHELP(W2,W,DT,FM1,2.,PS,PRAIN,E,RHOW,HF,Q,QM,...
            QE,HNEU);
    end
    UB2 = UB + DT/2.0*(Q1 + Q);
    niter=niter+1;
    if(niter >= imax)  %had * steps to converge now hit it.
        % What follows is a fix to numerical instability that results from
        % ninlinearity when the snowpack is shallow and melting rapidly.  If
        % convergence does not occur when the snowpack is not melting (a very
        % rare thing I just accept the predictor corrector solution.
        
        % The fix is asume the liquid fraction of the snow remains constant.
        % This to some extent bypasses the melt outflow estimates.
        % ae is added energy during the time step.
        AE=(Q1+Q+QM1+QM)*0.5*DT;
        %  This fix is only physically sensible under melting conditions
        %  and when ae is constant
        if((UB > 0) && (AE > 0) && (W > 0))
            %  Determine the w/ub ratio at the beginning of the time step.
            %  If liquid fraction is constant this is constant.
            %  Solve for ub2 and w2 that satisfy r=w2/ub2 and the added energy
            %  expression above.
            E2=(E+E1)*0.5;   % This is the average sublimation
            rlf=(UB+AE)/(RHOW*W*HF);
            if rlf>=1
                MR=W/DT+(PS+PRAIN-E2);
                if MR<0
                    MR=0;
                end
                QM=MR*RHOW*HF;
                W2=0;
                UB2=UB+AE-QM*DT;
            else
                r=W/UB;
                UB2=(RHOW*HF*(W+(PS+PRAIN-E2)*DT)-AE-UB)/(RHOW*HF*r-1);
                W2=r*UB2;
                %  Now changes in W imply changes in melt rate and melt advected
                %  energy.  These need to be adjusted so that mass and energy
                %  balances are still consistent.
                if W2<0
                    W2=0;
                    MR=(W-W2)/DT-E2+PS+PRAIN;
                end
                if MR<0
                    MR=0;
                    W2=W+(PS+PRAIN-E2)/DT;
                    if W2<0
                        W2=0;
                    end
                end
                QM=(AE+UB-UB2)/DT;
                MR=QM/(RHOW*HF);
                %                         QM=MR*RHOW*HF;
                %                         UB2=UB+AE-QM*DT;
            end
            
            Q=AE/DT-QM;
            %  Now fake the averages below
            QM1=QM;
            MR1=MR;
            Q1=Q;
        end
    end
    
end
W = W2;
UB=UB2;
% average values from two time steps for output.  This is done for mr
% and e to ensure mass balance and the others for better physical
% comparisons
QH=(QH+QH1)*0.5;
QE=(QE+QE1)*0.5;
E=(E+E1)*0.5;
MR=(MR+MR1)*0.5;
QM=(QM+QM1)*0.5;
TSURF=(TSURF+TSURF1)*0.5;
QNET=(QNET+QNET1)*0.5;
Q=(Q+Q1)*0.5;

end