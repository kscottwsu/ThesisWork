h = enthalpy(298);

hrxn = h.H2O - h.H2 - h.O2/2;

MHyd = 2/1000;
OMHyd = 0.3*MHyd;
Puse = -hrxn/0.7;
Puse = Puse*1e-6/(60*60);
LCOE = OMHyd/Puse;