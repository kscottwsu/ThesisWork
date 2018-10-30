function [LCOE,Penetration,iWind,iSolar,CurtailmentRate] = CurtailedRenewables(PercSolar,Scale,A)

[RenewablesData,CurtailableDemandData,iSolar,iWind,RampLimit] = ProfileGeneration(PercSolar,A);
MaxRamp = 5*RampLimit;%maxramp is in %/min, but 5 minute long intervals


iSolar = iSolar*Scale;
iWind = iWind*Scale;

RenewablesProfile = RenewablesData*Scale;
Excess = RenewablesProfile - CurtailableDemandData;
Excess(Excess<0) = 0;
UsedRenewables = RenewablesProfile - Excess;
Demand = CurtailableDemandData - RenewablesProfile;
Demand(Demand<0) = 0;

%ramp limitations accounted for with curtailment of renewables
% for i = 1:(length(Demand)-1)% demand data has irregularities, so there are a couple outliers
%     Ramp = Demand(i+1) - Demand(i);
%     if Ramp < -MaxRamp
%         Error = Ramp + MaxRamp;
%         Fix = max(-UsedRenewables(i+1),Error);
%         UsedRenewables(i+1) = UsedRenewables(i+1) + Fix;
%         Excess(i+1) = Excess(i+1) - Fix;
%         Demand(i+1) = Demand(i+1) - Fix;
%     end
%     if Ramp > MaxRamp
%         Error = Ramp - MaxRamp;
%         Fix = min(Excess(i+1),Error);
%         UsedRenewables(i+1) = UsedRenewables(i+1) + Fix;
%         Excess(i+1) = Excess(i+1) - Fix;
%         Demand(i+1) = Demand(i+1) - Fix;
%     end
% end

RenewablesDemandMet = sum(UsedRenewables);
UnmetDemand = sum(Demand);
CurtailedRenewables = sum(Excess);
CurtailmentRate = CurtailedRenewables/(CurtailedRenewables+RenewablesDemandMet);
Penetration = RenewablesDemandMet/(sum(CurtailableDemandData)/(1-A));


BaseLCOE = (iSolar*58.1*(0.27/0.18) + 44.3*iWind*(0.367/0.32))/(iWind+iSolar);
LCOE = BaseLCOE/(1-CurtailmentRate);
end

