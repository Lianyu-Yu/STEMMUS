%************************** hourlyRI () **********************
%               To get hourly radiation index

function [HRI,COSZEN]=hyri(MONTH,DAY,HOUR,DT,SLOPE,AZI,LAT)
    CRAD=pi()/180.0;
    %    CONVERT TIMES TO RADIANS FROM NOON
    T=(HOUR-12.0)*pi()/12.0;
    DELT1=DT*pi()/12.0;
    %    CONVERT TO RADIANS
    SLOPE1=SLOPE*CRAD;
    AZI1=AZI*CRAD;
    LAT1=LAT*CRAD;
    FJULIAN=JULIAN(MONTH,DAY);
    D=CRAD*23.5*sin((FJULIAN-82.0)*0.0172142);
    LP=asin(sin(SLOPE1)*cos(AZI1)*cos(LAT1)...
        +cos(SLOPE1)*sin(LAT1));
    TD=acos(-tan(LAT1)*tan(D));
    TPI=-tan(LP)*tan(D);
        if (abs(TPI)< 1.0)
            TP=acos(TPI);
        else
            TP=0.0;
        end
    DDT=atan(sin(AZI1)*sin(SLOPE1)/(cos(SLOPE1)*cos(LAT1)...
        -cos(AZI1)*sin(SLOPE1)*sin(LAT1)));
    T1=max(T,max(-TP-DDT,-TD));
    T2=min(T+DELT1,min(TD,TP-DDT));
    %     write(6,*)t1,t2
        if (T2<=T1)
            HRI=0.0;
        else
            HRI=(sin(D)*sin(LP)*(T2-T1)+cos(D)*cos(LP)*(sin(T2+DDT)...
                -sin(T1+DDT)))/(cos(SLOPE1)*DELT1);
        end
    COSZEN = HRI*cos(SLOPE1);
end