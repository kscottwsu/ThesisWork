function [LCOE] = LCOEcalc(iDemand,OandM,n,InitialCost,DegredationRate)
DiscountRate = 0.056;

TOandM = iDemand*OandM;

for i = 1:n
    if i == 1
        LCOE1 = InitialCost;
        LCOE2 = 0;
    end
    LCOE1 = LCOE1 + TOandM/(1+DiscountRate)^i;
    LCOE2 = LCOE2 + (iDemand*(1-DegredationRate)^i)/((1+DiscountRate)^i);
end
LCOE = LCOE1/LCOE2;
end