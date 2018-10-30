function [MRR] = ThermalLimitCalculation( )
load('BPA Thermal Data.mat');

ThermalMax = max(BPAData);
Ramp = BPAData(2:end)-BPAData(1:(end-1));
RampRate = (Ramp/ThermalMax)/5;
MinRampRate = min(RampRate);
MaxRampRate = max(RampRate);

MRR = max(-MinRampRate,MaxRampRate);
end

