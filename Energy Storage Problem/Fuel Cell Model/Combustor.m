function [ Inlet,OxyAdd ] = Oxidizer( Inlet )
OxyAdd.T = Inlet.T;
OxyAdd.O2 = Inlet.H2/2 + 2*Inlet.CH4;

RH2Combustion = Inlet.H2;
RCH4Combustion = Inlet.CH4;

h = enthalpy(298);
hrxn1 = h.H2O-h.H2-h.O2/2; %H2 + O2/2 -->  H2O, Hydrogen Combustion
hrxn2 = 2*h.H2O + h.CO2 -h.CH4 - 2*h.O2;% CH4 Combustion

delH = (-hrxn1*RH2Combustion - hrxn2*RCH4Combustion);
Inlet.H2 = 0;
Inlet.CH4 = 0;
Inlet.CO2 = Inlet.CO2 + RCH4Combustion;
Inlet.H2O = Inlet.H2O + 2*RCH4Combustion + RH2Combustion;

Cp = SpecHeat(Inlet);
Inlet.T = Inlet.T + delH/(Cp*NetFlow(Inlet));
end

