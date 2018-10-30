function [ Flow, QMethaneReform ] = Reformer( Flow )

h = enthalpy(298);
hrxn3 = 3*h.H2+h.CO-h.CH4-h.H2O; %CH4+H2O --> CO + 3H2, Methane reforming

rReform = 0.25;
RMethaneReform = rReform*Flow.CH4;
Flow.CH4 = Flow.CH4 - RMethaneReform;
if isfield(Flow, 'H2')
    Flow.H2 = Flow.H2 + 3*RMethaneReform;
else
    Flow.H2 = 3*RMethaneReform;
end
Flow.H2O = Flow.H2O - RMethaneReform;
if isfield(Flow, 'CO')
    Flow.CO = Flow.CO + RMethaneReform;
else
    Flow.CO = RMethaneReform;
end
QMethaneReform = RMethaneReform*hrxn3/1000;%heat demand for external methane reformer
end

