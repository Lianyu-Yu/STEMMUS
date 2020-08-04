%Computes the surface temperature of snow
function [SRFTMPO]=SRFTMPO(QSI,A,QLI,QPIN,EA,TA,TAVE,TK,PR,RA,CP,RHO,RKN,...
    HNEU,ES,SBC,CS,RS,W,qnetob,IRADFL,WS,Z,G,FC,fstab)
%     dimension tint(0:10,2)
    ncall=0;tol=0.05;
    ncall = ncall+1;
    NITER = 10;
    TSTAR = TA;
    QP=QPIN;    %store input variable locally without changing global value
    if(W<=0 && QP>0),QP=0;end
    %      if(w.le.0. and. qp.gt.0.)qp=0.   %ignore the effect of precip advected
    %     energy on the calculation of surface temperature when there is no snow.
    %     Without this ridiculously high temperatures can result as the model
    %     tries to balance outgoing radiation with precip advected energy.
    TSURF = SNOTMP(TSTAR,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,...
        RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab);
    %     The first calculation gets the solution on the right side of ta
    %     since it uses neutral stability and the linearization (equivalent
    %     to Newton-Rhapson) will move in the direction of a solution in the case
    %     of a well behaved function.  See Notebook 8/3/94 for further elaboration.
    %     Use this to place bounds on surface temperature for iteration.
    if TSURF>TA
        tlb=TA;
        tub=TA+30;   % Max upper bound 30 C warmer than surface
        if(TSURF>tub),TSURF=(tlb+tub)*0.5;end  %sometimes TSURF is outside these bounds
    else
        tlb=TA-30;
        tub=TA;
        if(TSURF<tlb),TSURF=(tlb+tub)*0.5;end
    end
    tint(1,1)=TSTAR;
    tint(1,2)=TSURF;
    %  Now iterate
    TSTAR=TSURF;
    for i=1:NITER
        TSURF = SNOTMP(TSTAR,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,...
            RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab);
        
        tint(i,1)=TSTAR;
        tint(i,2)=TSURF;
        if (tlb<=TSURF && TSURF <= tub)
            if(TSURF>TSTAR)
                tlb=TSTAR;   %increasing so can increase lower bound.
            else
                tub=TSTAR;
            end
        elseif(TSURF > tub)   %upper bound overshot
            tlb=TSTAR;                   %increase lower bound and
            TSURF=(TSTAR+tub)*0.5;        %guess at solution halfway
        else    %TSURF .lt. tlb  here, i.e. lower bound overshot
            tub=TSTAR;
            TSURF=(TSTAR+tlb)*0.5;
        end
        %   Check for convergence
        if(abs(TSURF-TSTAR)<tol)
            break; %GO TO 20
        else
            TSTAR = TSURF;
            %   Newton rhapson not converging so use bisection
            f1=sntmpb(tlb,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,...
                RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab);
            f2=sntmpb(tub,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,...
                RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab);
            TSURF=(tlb+tub)*0.5;
            if f1*f2>0
                print 'SRFTMP has failed to find a solution' %,ncall...
                % ,tlb,tub,TSURF
            else
                nib=(log(tub-tlb)-log(tol))/log(2);
                for j=1:nib
                    f=sntmpb(TSURF,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,...
                        RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab);
                    if f*f1>0
                        tlb=TSURF;
                        f1=f;
                    else
                        tub=TSURF;
                        %         f2=f
                    end
                    TSURF=(tlb+tub)*0.5;
                end
            end
            %         end
        end
       
        if (W>0 && TSURF>0)
            SRFTMPO = 0;
        else
            SRFTMPO = TSURF;
        end 
    end
end