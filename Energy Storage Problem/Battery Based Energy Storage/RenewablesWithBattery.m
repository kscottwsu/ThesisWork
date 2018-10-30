function [LCOETotal,Penetration,iWind,iSolar] = RenewablesWithBattery(PercSolar,Scale,A,BatterySizing)
RTE = 0.81;%round trip efficiency
BinEff = RTE^0.5;
BoutEff = RTE^0.5;
BatteryDecay = 0.9995434;

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
CurtailedRenewables = sum(Excess);
CurtailmentRate = CurtailedRenewables/(CurtailedRenewables+RenewablesDemandMet);

BaseLCOE = (iSolar*58.1*(0.27/0.18)+44.3*iWind*(0.367/0.32))/(iWind+iSolar);
LCOEnoStorage = BaseLCOE/(1-CurtailmentRate);

%SmoothedDemand = Smoothing(Demand,144,1);

DemandBattery = Demand ;%- SmoothedDemand;%find out where you can peak shave
DemandBattery(DemandBattery<0) = 0;
DemandBattery2 = DemandBattery;
Excess2 = Excess;

%DemandBattery = Demand;
BatteryProfile = 0*Demand;
BatteryStorage = zeros(length(Demand)+1,1);


MaxCRate = 1;
MaxRate = MaxCRate*BatterySizing*5/60;


% Battery calculation
for k = 1:3
    DemandBattery = DemandBattery2;
    Excess = Excess2;
    if k == 1
        BatteryStorage(1) = 0;
    else
        BatteryStorage(1) = BatteryStorage(end);
    end
    for i = 1:(length(Demand)-1)
        BatteryStorage(i+1) = BatteryStorage(i)*BatteryDecay;
        if Excess(i) > 0%doesn't check regions where renewables were alread curtailed
            temp = min([Excess(i),(BatterySizing-BatteryStorage(i+1))/BinEff,MaxRate/BinEff]);
            BatteryUse = temp*BinEff;
            BatteryStorage(i+1) = BatteryStorage(i+1) + BatteryUse;
            BatteryProfile(i) = -temp;
            Excess(i) = Excess(i) - temp;
        else
            BatteryUse = min([DemandBattery(i)/BoutEff,BatteryStorage(i+1),MaxRate/BoutEff]);%this is where forecasting would go
            BatteryStorage(i+1) = BatteryStorage(i+1) - BatteryUse;
            DemandBattery(i) = DemandBattery(i) - BatteryUse*BoutEff;
            BatteryProfile(i) = BatteryUse;
        end
    end
end

BatterySupplied = DemandBattery2 - DemandBattery;
BatteryDemandMet = sum(BatterySupplied);
Demand = Demand - BatterySupplied;

cycles = BatteryDemandMet/BatterySizing;
if BatteryDemandMet == 0 && BatterySizing == 0
    LCOETotal = LCOEnoStorage;
elseif BatteryDemandMet == 0
    LCOETotal = LCOEnoStorage + BatterySizing*250*1000/RenewablesDemandMet;
else
    LCOEBattery = LCOEcalc(BatteryDemandMet,2,15,BatterySizing*250*1000,0.01);%based on eight hour utility battery storage capictal costs, nrel
    LCOETotal = (LCOEBattery*BatteryDemandMet + LCOEnoStorage*RenewablesDemandMet)/(RenewablesDemandMet+BatteryDemandMet);
end

CurtailedRenewables = sum(Excess);
UnMetDemand = sum(Demand);
MetDemand = sum(CurtailableDemandData -Demand);

CurtailmentRate = CurtailedRenewables/(iWind+iSolar);
Penetration = MetDemand/(sum(CurtailableDemandData)/(1-A));
NoStoragePenetraion = RenewablesDemandMet/(sum(CurtailableDemandData)/(1-A));
end