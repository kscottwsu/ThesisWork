function [ LCOE,Penetration,iWind,iSolar,RoundTripEfficiency,StorageSize,MasHyd,CurtailmentRate,NoStoragePenetration,MasHHTM,MoneyHyd ] = RenewablesWithFuelCell( PercSolar,Scale,A,FCSizing )

[RenewablesData,CurtailableDemandData,iSolar,iWind,RampLimit] = ProfileGeneration(PercSolar,A);
MaxRamp = 5*RampLimit;

iSolar = iSolar*Scale;
iWind = iWind*Scale;

RenewablesProfile = RenewablesData*Scale;
Excess = RenewablesProfile - CurtailableDemandData;
Excess(Excess<0) = 0;
UsedRenewables = RenewablesProfile - Excess;
Demand = CurtailableDemandData - RenewablesProfile;
Demand(Demand<0) = 0;

RenewablesDemandMet = sum(UsedRenewables);
UnmetDemand = sum(Demand);
CurtailedRenewablesNoStorage = sum(Excess);
CurtailmentRateNoStorage = CurtailedRenewablesNoStorage/(CurtailedRenewablesNoStorage+RenewablesDemandMet);

BaseLCOE = (iSolar*58.1*(0.27/0.18)+44.3*iWind*(0.367/0.32))/(iWind+iSolar);
LCOEnoStorage = BaseLCOE/(1-CurtailmentRateNoStorage);

DemandFC = Demand ;%- SmoothedDemand;%find out where you can peak shave
DemandFC(DemandFC<0) = 0;
DemandFC2 = DemandFC;
Excess2 = Excess;

FCProfile = 0*Demand;
FCStorage = zeros(length(Demand)+2,1);

load('EffCurve.mat');
EffCurve.MaxFCPower = 0.515;
EffCurve.MinFCPower = -1.3;

eFCSizing = 0.99*FCSizing;%assume 99% uptime


FCElectrolysisMax = -eFCSizing*EffCurve.MinFCPower;
FCMax = eFCSizing*EffCurve.MaxFCPower;


h = enthalpy(298);
HMethane = -(2*h.H2O + h.CO2 - h.CH4 - 2*h.O2);
HHyd = h.H2 + 0.5*h.O2 - h.H2O;



Ein = 0;
Eout = 0;
for k = 1:2
    DemandFC = DemandFC2;
    Excess = Excess2;
    HydGeneration = 0*Demand;
    HTM = 0;
    Ein = 0;
    Eout = 0;
    
    if k == 1
        FCStorage(1) = 0;
    else
        FCStorage(1) = FCStorage(end);
    end
    for i = 1:(length(Demand))
        FCStorage(i+1) = FCStorage(i);%stored energy losses applied here
        if Excess(i) > 0
            Pin = min([Excess(i)*(60*60)/300,FCElectrolysisMax])/eFCSizing;%Power Density in Fuel Cell
            [Eff, ~] = Efficiency(-Pin,EffCurve);
            PStored = Eff*Pin*eFCSizing;
            MolCH4 = 1e6*PStored*300/HMethane;
            FCStorage(i+1) = FCStorage(i+1) + MolCH4;%in  MolCH4
            FCProfile(i) = -Pin*eFCSizing;%in MW
            Excess(i) = Excess(i) - (Pin*eFCSizing)*300/(60*60);%in MWh
            Ein = Ein + 300*eFCSizing*Pin/(60*60);%in MWh
        else
            EStorage = FCStorage(i+1)*HMethane*1e-6;%in MJ
            Pout = min([DemandFC(i)*(60*60)/300,FCMax])/eFCSizing;%power density in MW
            [Eff, HydProduction] = Efficiency(Pout,EffCurve);
            if EStorage > 0
                PStored = min(Pout*eFCSizing/Eff,EStorage/300);%MW
                if EStorage/300 < Pout*eFCSizing/Eff
                    HydProduction = HydProduction*(EStorage/300)/(Pout*eFCSizing/Eff);
                end
            else
                PStored = 0;
            end
            MolCH4 = 300*PStored/(HMethane*1e-6);
            FCStorage(i+1) = FCStorage(i+1) - MolCH4;
            FCProfile(i) = PStored*Eff;%MW
            DemandFC(i) = DemandFC(i) - PStored*Eff*300/(60*60);%mWh
            Eout = Eout + 300*PStored*Eff/(60*60);%in MWh
            if HydProduction < 0
                HydProduction = 0;
            end
            HydGeneration(i) = 300*1e6*eFCSizing*HydProduction;%in Mols H2
            HTM = HTM + HydGeneration(i);
            HydProduction = 0;
        end
    end
    FCStorage(end) = FCStorage(end-1);
    FCStorage = FCStorage - min(FCStorage);
end

%if it gains 
if FCStorage(end)  > (FCStorage(1)+1)
    FCStorageMax = FCStorage(1) + max(FCStorage) - FCStorage(end);

    for k = 1:2
        DemandFC = DemandFC2;
        Excess = Excess2;
        HydGeneration = 0*Demand;
        HTM = 0;
        Ein = 0;
        Eout = 0;
        if k == 1
            FCStorage(1) = 0;
        else
            FCStorage(1) = FCStorage(end);
        end
        for i = 1:(length(Demand))
            FCStorage(i+1) = FCStorage(i);%stored energy losses applied here
            if Excess(i) > 0  && FCStorage(i+1) ~= FCStorageMax
                MaxMolCH4 = FCStorageMax - FCStorage(i+1);
                               
                Pin = min([Excess(i)*(60*60)/300,FCElectrolysisMax])/eFCSizing;%Power Density in Fuel Cell
                [Eff, ~] = Efficiency(-Pin,EffCurve);
                PStored = Eff*Pin*eFCSizing;
                MolCH4 = 1e6*PStored*300/HMethane;
                
                MolCH4 = min([MolCH4, MaxMolCH4]);
                P = MolCH4*1e-6*HMethane/300;%Power in MW
                Pin = P/(Eff*eFCSizing);
                
                FCStorage(i+1) = FCStorage(i+1) + MolCH4;%in  MolCH4
                FCProfile(i) = -Pin*eFCSizing;%in MW
                Excess(i) = Excess(i) - (Pin*eFCSizing)*300/(60*60);%in MWh
                Ein = Ein + 300*eFCSizing*Pin/(60*60);%in MWh
            elseif Excess(i) > 0  && FCStorage(i+1) == FCStorageMax
                Pin = min([Excess(i)*(60*60)/300,FCElectrolysisMax])/eFCSizing;
                Excess(i) = Excess(i) - Pin*eFCSizing*300/(60*60);
                [~, HydProduction] = EfficiencyHydElectrolysis(-Pin,EffCurve);
                if HydProduction < 0
                    HydProduction = 0;
                end
                HydGeneration(i) = 300*1e6*eFCSizing*HydProduction;
                HydProduction = 0;
            else
                EStorage = FCStorage(i+1)*HMethane*1e-6;%in MJ
                Pout = min([DemandFC(i)*(60*60)/300,FCMax])/eFCSizing;%power density in MW
                [Eff, HydProduction] = Efficiency(Pout,EffCurve);
                if EStorage > 0
                    PStored = min(Pout*eFCSizing/Eff,EStorage/300);%MW
                    if EStorage/300 < Pout*eFCSizing/Eff
                        HydProduction = 0*HydProduction*(EStorage/300)/(Pout*eFCSizing/Eff);
                    end
                else
                    PStored = 0;
                end
                MolCH4 = 300*PStored/(HMethane*1e-6);
                FCStorage(i+1) = FCStorage(i+1) - MolCH4;
                FCProfile(i) = PStored*Eff;%MW
                DemandFC(i) = DemandFC(i) - PStored*Eff*300/(60*60);%mWh
                Eout = Eout + 300*PStored*Eff/(60*60);%in MWh
                if HydProduction < 0
                    HydProduction = 0;
                end
                HydGeneration(i) = 300*1e6*eFCSizing*HydProduction;%in Mols H2
                HTM = HTM + HydGeneration(i);
                HydProduction = 0;
            end
        end
        FCStorage(end) = FCStorage(end-1);
        FCStorage = FCStorage - min(FCStorage);
    end


end

PUseEC = sum(abs(FCProfile(FCProfile<0)))*300/(60*60);%in MWh




StorageSize = max(FCStorage);
if isnan(StorageSize)
    StorageSize = 0;
end

if any(isnan(DemandFC))
    DemandFC = DemandFC2;
end


FCSupplied = DemandFC2 - DemandFC;
FCDemandMet = sum(FCSupplied);
Demand = Demand - FCSupplied;


%%

molCH4 = max(FCStorage);
R = 8.314;
T = 333.15;
p = 160*100*1000;%160 bar in pascals
VolCH4 = molCH4*R*T/p;
CH4StorageCost = 36e6*1.5*VolCH4/(800e3);
CO2StorageCost = CH4StorageCost;
O2StorageCost = 2*CH4StorageCost;

molH2O = 10*molCH4;
massH2O = molH2O*18/1000;%18 g/mol
volH2O = massH2O/1000;
WaterStorageCost = volH2O*50.4;

StorageCost = CH4StorageCost + CO2StorageCost + 4*O2StorageCost + WaterStorageCost;
LCOEStorage = LCOEcalc(FCDemandMet,1,20,StorageCost,0);

%safety factor = 50%
FCPower = 2*EffCurve.MaxFCPower;

%% based on wendel and jensen paper
FuelCellCost = FCPower*150*1000*FCSizing;%DOE Target
BOPCost = FCPower*(1075-200 - 155 - 128 - 23)*1000*FCSizing;%takes jensen's results and subtracts parts accounted for elsewhere

MasHHTM = HTM*2/1000;
CostHTM = MasHHTM*2.00;
effectiveLCOEHTM = CostHTM/FCDemandMet;%will generate hydrogen for 10 years

%hydrogen from electrolysis
MasHyd = sum(HydGeneration)*2/1000 - MasHHTM;%2 grams per mol
CostHyd = MasHyd*2.00;%$2 per gallon of gas equivalent, 1kg H2 = 1 gallon of gas equivalent
effectiveLCOEHyd = CostHyd/FCDemandMet;%will generate hydrogen for 10 years


MoneyHyd = CostHyd + CostHTM;

if MasHyd < 0
    MasHyd = 0;
end

if FCSizing == 0
    LCOETotal = LCOEnoStorage;
    LCOEHyd = 0;
    MasHyd = 0;
    MasHHTM = 0;
elseif FCDemandMet == 0
    MasHyd = 0;
    MasHHTM = 0;
    LCOEHyd = 0;
    LCOETotal = LCOEnoStorage + FCSizing*1075*1000/RenewablesDemandMet;
else
    LCOEHTM = effectiveLCOEHTM;
    LCOEHyd = effectiveLCOEHyd+LCOEHTM;
    LCOEOMFC = LCOEcalc(FCDemandMet,20,9,0,0.0263);
    LCOEOMEC = LCOEcalc(PUseEC,10,9,0,0.0263);
    LCOEFC = LCOEcalc(FCDemandMet,0,9,FuelCellCost,0.0263);%degredation/lifetime based on DOE 2020 target
    LCOEBOP = LCOEcalc(FCDemandMet,1,20,BOPCost,0.01);
    LCOETotal = (LCOEFC*FCDemandMet + LCOEBOP*FCDemandMet + LCOEOMFC*FCDemandMet + LCOEOMEC*PUseEC + LCOEStorage*FCDemandMet - LCOEHyd*FCDemandMet + LCOEnoStorage*RenewablesDemandMet)/(RenewablesDemandMet + FCDemandMet);
end


CurtailedRenewables = sum(Excess);
UnMetDemand = sum(Demand);
MetDemand = sum(CurtailableDemandData - Demand);

CurtailmentRate = CurtailedRenewables/(iWind+iSolar);
Penetration = MetDemand/(sum(CurtailableDemandData)/(1-A));
NoStoragePenetration = RenewablesDemandMet/(sum(CurtailableDemandData)/(1-A));
RoundTripEfficiency = Eout/Ein;
if isnan(Eout)
    RoundTripEfficiency = 0;
elseif FCSizing == 0
    RoundTripEfficiency = 0;
elseif isnan(RoundTripEfficiency)
    RoundTripEfficiency = 0;
end
LCOE = LCOETotal;
end

function [Eff, HydProduction] = Efficiency(P,EffCurve)

n1 = find(EffCurve.Power > P,1,'First');
if isempty(n1)
    n1 = length(EffCurve.Power);
end
if n1 == 1
    n1 = 2;
end
n2 = n1-1;
    
Eff = EffCurve.Eff(n2) +(EffCurve.Eff(n1)-EffCurve.Eff(n2))*(P-EffCurve.Power(n2))/(EffCurve.Power(n1)-EffCurve.Power(n2));
HydProduction = EffCurve.HydProduction(n2) +(EffCurve.HydProduction(n1)-EffCurve.HydProduction(n2))*(P-EffCurve.Power(n2))/(EffCurve.Power(n1)-EffCurve.Power(n2));
end
    
    
function [Eff, HydProduction] = EfficiencyHydElectrolysis(P,EffCurve)
n1 = find(EffCurve.ElectrolysisModePower < P,1,'First');
if isempty(n1)
    n1 = length(EffCurve.ElectrolysisModePower);
elseif n1 == 1
    n1 = 2;
end
n2 = n1-1;
HydProduction = EffCurve.ElectrolysisModeHyd(n2) + (EffCurve.ElectrolysisModeHyd(n1)-EffCurve.ElectrolysisModeHyd(n2))*(P-EffCurve.ElectrolysisModePower(n2))/(EffCurve.ElectrolysisModePower(n1)-EffCurve.ElectrolysisModePower(n2));
Eff = EffCurve.ElectrolysisModeEff(n1) + (EffCurve.ElectrolysisModeEff(n2)-EffCurve.ElectrolysisModeEff(n1))*(P-EffCurve.ElectrolysisModePower(n1))/(EffCurve.ElectrolysisModePower(n2)-EffCurve.ElectrolysisModePower(n1));
end