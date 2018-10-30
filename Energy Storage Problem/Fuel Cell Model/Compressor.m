function [ Outlet, Wc ] = Compressor( Inlet, Pin, Pout)
eff = 0.8;
Ru = 8.3145;
Outlet = Inlet;
Cp = SpecHeat(Outlet);
Gamma = Cp/(Cp-Ru);
T2s = Inlet.T*(Pout/Pin)^((Gamma -1)/Gamma);
[~,H1] = enthalpy(Inlet);%initial enthalpy
Outlet.T = T2s;
[~,H2s] = enthalpy(Outlet);
H2a = H1+(H2s-H1)/eff;%the real enthalpy after compression
Wc = (H2a - H1)/1000;%in KW
dHa = (H2a - H2s);%in W
Cp = SpecHeat(Outlet);
dTa = dHa/(Cp*NetFlow(Inlet));%assume all entropy increase goes into gas
Outlet.T = Outlet.T + dTa;
end

