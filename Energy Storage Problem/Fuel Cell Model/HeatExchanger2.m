function [ InletFlow1, InletFlow2, QHX ] = HeatExchanger2( InletFlow1,InletFlow2 )
%assume 100% heat exchanger efficiency

T1 = InletFlow1.T;
T2 = InletFlow2.T;

delT = abs(T1 - T2);
InletFlow1.T = (T1+T2)/2;
InletFlow2.T = (T1+T2)/2;


Cp1 = SpecHeat(InletFlow1);
Cp2 = SpecHeat(InletFlow2);
N1 = NetFlow(InletFlow1);
N2 = NetFlow(InletFlow2);

EMax1 = delT*Cp1*N1;
EMax2 = delT*Cp2*N2;

QHX = min(EMax1,EMax2);

if EMax2 > EMax1
    InletFlow1.T = T2;
    InletFlow2.T = T2 + Cp1*N1*(T1-T2)/(Cp2*N2);
elseif EMax1 > EMax2
    InletFlow1.T = T1 + Cp2*N2*(T2-T1)/(Cp1*N1);
    InletFlow2.T = T1;
elseif EMax1 == EMax2
    InletFlow1.T = T2;
    InletFlow2.T = T1;
end

end