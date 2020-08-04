%C***************************** JULIAN () ****************************
%             To convert the real date to julian date

function [JULIAN]=JULIAN(MONTH,DAY)
      MADD(1:12)=[0,31,59,90,120,151,181,212,243,273,304,334];
      JULIAN=DAY+MADD(MONTH);
end
