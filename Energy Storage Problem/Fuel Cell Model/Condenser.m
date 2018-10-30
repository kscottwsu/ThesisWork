function [ Flow, Q ] = Condenser( Flow, T )
Eff = 0.9;
Hv = 40.660;
Qv = Flow.H2O*Hv;

[~,H1] = enthalpy(Flow);
Flow.T = T;
[~,H2] = enthalpy(Flow);


Q = Qv + (H2-H1)/1000;
Q = Q*Eff;

Flow.H2O = 0;
end

