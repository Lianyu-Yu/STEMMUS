%****************************** QLIF () *******************************
%    Computes the incoming longwave radiation using satterlund Formula
%     Modified 10/13/94 to account for cloudiness.
%     Emissivity of cloud cover fraction is assumed to be 1.
%
function [QLIFF]=qlif(TA,RH,TK,SBC,cf)
      EA = SVPW(TA)*RH;
      TAK = TA + TK;
      EA1 = 1.08*(1.0-exp(-(EA/100.0)^((TAK)/2016.0)));
      QLIFF =(cf+(1-cf)*EA1)*SBC*TAK^4;
end