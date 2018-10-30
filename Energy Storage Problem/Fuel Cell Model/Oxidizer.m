function [ Inlet,OxyAdd ] = Oxidizer( Inlet,TMax )
OxyAdd.T = 333.15;
OxyAdd.O2 = Inlet.H2/2 + 2*Inlet.CH4 + Inlet.CO;

RH2Combustion = Inlet.H2;
RCH4Combustion = Inlet.CH4;
RCOCombustion = Inlet.CO;

h = enthalpy(298);
hrxn1 = h.H2O-h.H2-h.O2/2; %H2 + O2/2 -->  H2O, Hydrogen Combustion
hrxn2 = 2*h.H2O + h.CO2 -h.CH4 - 2*h.O2;% CH4 Combustion
%hrxn3 = h.CO2 - h.CO - 0.5*h.O2;

delHOxy = enthalpy(OxyAdd);
OxyAdd.T = Inlet.T;
delHOxy = delHOxy - enthalpy(OxyAdd);
delH = (-hrxn1*RH2Combustion - hrxn2*RCH4Combustion) + delHOxy;

Inlet.H2 = 0;
Inlet.CH4 = 0;
Inlet.CO = 0;
Inlet.CO2 = Inlet.CO2 + RCH4Combustion + RCOCombustion;
Inlet.H2O = Inlet.H2O + 2*RCH4Combustion + RH2Combustion;

Cp = SpecHeat(Inlet);
Inlet.T = Inlet.T + delH/(Cp*NetFlow(Inlet));

if Inlet.T > TMax
    Inlet.T = TMax;
end
end

