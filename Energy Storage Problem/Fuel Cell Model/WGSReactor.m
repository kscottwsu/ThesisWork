function [ Flow, QWGS] = WGSReactor( Flow )
h = enthalpy(298);
hrxn = h.CO2+h.H2-h.CO-h.H2O; %CO + H20 --> CO2 + H2, Water Gas Shift

if ~isfield(Flow,'CO')
    Flow.CO = 0;
end
if ~isfield(Flow,'CO2')
    Flow.CO2 = 0;
end
if ~isfield(Flow,'CH4')
    Flow.CH4 = 0;
end
if ~isfield(Flow,'H2O')
    Flow.H2O = 0;
end


rWGS = 1;
RWGS = rWGS*Flow.CO;
RWGS = min([RWGS,Flow.CO,Flow.H2O]);
RWGS = max(RWGS,0);
Flow.CO = Flow.CO - RWGS;
Flow.H2O = Flow.H2O - RWGS;
if isfield(Flow,'CO2')
    Flow.CO2 = Flow.CO2 + RWGS;
else
    Flow.CO2 = RWGS;
end
Flow.H2 = Flow.H2 + RWGS;

QWGS = -RWGS*hrxn/1000;
end