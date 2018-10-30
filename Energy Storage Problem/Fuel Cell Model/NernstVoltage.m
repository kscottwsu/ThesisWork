function [ E0, V_Nernst ] = NernstVoltage(X,P)
F = 96485.332896;%faraday constant
Ru = 8.3145;%Ideal Gas Constant

T = X.T;

[h,~] = enthalpy(T);
s = entropy(T);
E0 = -((h.H2O-s.H2O.*T)-(h.H2-s.H2.*T)-.5*(h.O2-s.O2.*T))/(2*F);%ideal voltage
EThermoNeutral = -(h.H2O-h.H2-0.5*h.O2)/(2*F);
frac = sqrt(X.O2)*X.H2/X.H2O;
V_Nernst = Ru*T/(2*F)*log(abs(sqrt(P/101.325)*frac));%factors in effects of partial pressures

end

