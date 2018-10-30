function [ Flow,Q ] = ChemicalEquilibrium(Flow,P)
%P: pressure in bar

R = 8.3145;
T = Flow.T;
Flow = FlowInitialize(Flow);

h = enthalpy(298);
s = entropy(298);

hSR = 3*h.H2+h.CO-h.CH4-h.H2O;
hCD = 2*h.H2+h.C-h.CH4;
hWGS = h.H2 + h.CO2 - h.CO - h.H2O;

sSR = 3*s.H2+s.CO-s.CH4-s.H2O;
sCD = 2*s.H2+s.C-s.CH4;
sWGS = s.H2+s.CO2-s.CO-s.H2O;

gSR = hSR-sSR*T;
gCD = hCD - sCD*T;
gWGS = hWGS - sWGS*T;

KSR = exp(-gSR/(R*T));
KCD = exp(-gCD/(R*T));
KWGS = exp(-gWGS/(R*T));

x(1) = 0;%steam reforming reaction rate
x(2) = 0;%carbon deposition reaction rate
x(3) = 0;%water gas shift reaction rate
Temp = Flow;

%%
imax = 20;
for i = 1:imax
NF = NetFlow(Flow);
PP.H2 = P*Flow.H2/NF;
PP.H2O = P*Flow.H2O/NF;
PP.CH4 = P*Flow.CH4/NF;
PP.CO = P*Flow.CO/NF;
PP.CO2 = P*Flow.CO2/NF;
PP.C = P*Flow.C/NF;


dx = 1e-6*mean([Flow.H2,Flow.H2O,Flow.CH4,Flow.CO,Flow.CO2]);

Fp(1) = KSR*PP.CH4*PP.H2O - PP.CO*PP.H2^3;
Fp(2) = KCD*PP.CH4 - PP.C*PP.H2^2;
Fp(3) = KWGS*PP.CO*PP.H2O - PP.CO2*PP.H2;

%steam reforming
Flow2 = Flow;
Flow2.H2O = Flow2.H2O - dx;
Flow2.CH4 = Flow2.CH4 - dx;
Flow2.H2 = Flow2.H2 + 3*dx;
Flow2.CO = Flow2.CO + dx;

NF = NetFlow(Flow2);
PP.H2 = P*Flow2.H2/NF;
PP.H2O = P*Flow2.H2O/NF;
PP.CH4 = P*Flow2.CH4/NF;
PP.CO = P*Flow2.CO/NF;
PP.CO2 = P*Flow2.CO2/NF;
PP.C = P*Flow2.C/NF;

J(1,1) = ((KSR*PP.CH4*PP.H2O - PP.CO*PP.H2^3)-Fp(1))/dx;
J(2,1) = ((KCD*PP.CH4 - PP.C*PP.H2^2)-Fp(2))/dx;
J(3,1) = ((KWGS*PP.CO*PP.H2O - PP.CO2*PP.H2)-Fp(3))/dx;

%carbon deposition
Flow2 = Flow;
Flow2.CH4 = Flow2.CH4 - dx;
Flow2.H2 = Flow2.H2 + 2*dx;
Flow2.C = Flow2.C + dx;

NF = NetFlow(Flow2);
PP.H2 = P*Flow2.H2/NF;
PP.H2O = P*Flow2.H2O/NF;
PP.CH4 = P*Flow2.CH4/NF;
PP.CO = P*Flow2.CO/NF;
PP.CO2 = P*Flow2.CO2/NF;
PP.C = P*Flow2.C/NF;

J(1,2) = ((KSR*PP.CH4*PP.H2O - PP.CO*PP.H2^3)-Fp(1))/dx;
J(2,2) = ((KCD*PP.CH4 - PP.C*PP.H2^2)-Fp(2))/dx;
J(3,2) = ((KWGS*PP.CO*PP.H2O - PP.CO2*PP.H2)-Fp(3))/dx;

%WGS
Flow2 = Flow;
Flow2.CO = Flow2.CO - dx;
Flow2.H2O = Flow2.H2O - dx;
Flow2.H2 = Flow2.H2 + dx;
Flow2.CO2 = Flow2.CO2 + dx;

NF = NetFlow(Flow2);
PP.H2 = P*Flow2.H2/NF;
PP.H2O = P*Flow2.H2O/NF;
PP.CH4 = P*Flow2.CH4/NF;
PP.CO = P*Flow2.CO/NF;
PP.CO2 = P*Flow2.CO2/NF;
PP.C = P*Flow2.C/NF;


J(1,3) = ((KSR*PP.CH4*PP.H2O - PP.CO*PP.H2^3)-Fp(1))/dx;
J(2,3) = ((KCD*PP.CH4 - PP.C*PP.H2^2)-Fp(2))/dx;
J(3,3) = ((KWGS*PP.CO*PP.H2O - PP.CO2*PP.H2)-Fp(3))/dx;

x = x - (J^-1*Fp')';

Flow.CO = Temp.CO + x(1) + 2*x(2) - x(3);
Flow.H2O = Temp.H2O - x(1) - x(3);
Flow.H2 = Temp.H2 + 3*x(1) + 2*x(2) + x(3);
Flow.CO2 = Temp.CO2 + x(3);
Flow.CH4 = Temp.CH4 -x(1) - x(2);
Flow.C = Temp.C + x(2);
% end
% NF = NetFlow(Flow);
% MF.CO = Flow.CO/NF;
% MF.CO2 = Flow.CO2/NF;
% MF.CH4 = Flow.CH4/NF;
% MF.H2 = Flow.H2/NF;
% MF.H2O = Flow.H2O/NF;

Q = (-hSR*x(1) + -hCD*x(2) + -hWGS*x(3))/1000;
end

end
