 %Partitioning of precipitation into rain and snow
      
 function [PARTSNOW]=PARTSNOW(P,TA,TR,TS)

      if (TA<TS) 
        PARTSNOW=P;
      elseif (TA>TR) 
        PARTSNOW=0.0;
      else
        PARTSNOW=P*(TR-TA)/(TR-TS);
      end

 end