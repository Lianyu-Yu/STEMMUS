%C************************** UPDATETIME () ***************************
%                Update Time for each time step

function [YEAR,MONTH,DAY,HOUR]=UPDATETIME(YEAR,MONTH,DAY,HOUR,DT)
%       DIMENSION DMON(12)
%       INTEGER YEAR,DAY,DMON,DM
      DMON(1:12)=[31,28,31,30,31,30,31,31,30,31,30,31];
      HOUR=HOUR+DT;
      if (mod(YEAR,4)>0) 
         DMON(2)=28;
      else
         DMON(2)=29;
      end
      DM=DMON(MONTH);
      if (HOUR>=24.0) 
        HOUR=HOUR-24.0;
        DAY=DAY+1;
      else
        HOUR=HOUR;
      end
      if (DAY>DM) 
        DAY=1;
        MONTH=MONTH+1;
      else
        MONTH=MONTH;
      end
      if (MONTH>12) 
        MONTH=1;
        YEAR=YEAR+1;
      else
        YEAR=YEAR;
      end
      DMON(2)=28;
%       RETURN
end
