function [ Water, QEvaporator ] = Evaporator( Water )
Eff = 0.9;
%2257;%KJ/Kg
Hv = 40.660;%KJ/mol heat of vaporization
Qv = Water.H2O*Hv;

[~,H1] = enthalpy(Water);
Water.T = 373.15;%boiling point in K
[~,H2] = enthalpy(Water);

QEvaporator = Qv + (H2-H1)/1000;
QEvaporator = QEvaporator/Eff;
end

