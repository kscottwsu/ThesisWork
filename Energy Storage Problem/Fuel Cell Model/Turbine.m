function [ Outlet,Wt ] = Turbine( Inlet, Pin, Pout )
eff = 0.88;
Ru = 8.3145;
Outlet = Inlet;
Cp = SpecHeat(Outlet);
Gamma = Cp/(Cp-Ru);
T2s = Inlet.T*(Pout/Pin)^((Gamma -1)/Gamma);
[~,H1] = enthalpy(Inlet);%initial enthalpy
Outlet.T = T2s;
[~,H2s] = enthalpy(Outlet);
H2a = H1 - (H1-H2s)*eff;%the real enthalpy after Turbine
Wt = (H1-H2a)/1000;
dHa = (H2a - H2s);
dTa = dHa/(Cp*NetFlow(Inlet));%assume all work goes into gas
Outlet.T = Outlet.T + dTa;
end