function [Flow] = FlowInitialize(Flow)
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
if ~isfield(Flow,'H2')
    Flow.H2 = 0;
end
if ~isfield(Flow,'C')
    Flow.C = 0;
end
end