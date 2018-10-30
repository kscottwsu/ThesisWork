function [  ] = OptimalFuelCellPlotting(  )
%Calculates Optimal LCOE with Reversible SOFC based energy storage
A = 0.5;%percentage of demand met by baseline supply.
imax = 70;
jmax = 70;
maxScale = 1.4;
tol = 1e-5;

load('EffCurve.mat');
EffCurve.MaxFCPower = 0.515;
EffCurve.MinFCPower = -1.3;

for j = 1:jmax
    tic
    PercSolar(j) = 0.35*((j-1)/(jmax-1));
    for i = 1:imax
        Scale(i) = (maxScale)*(i/imax);
        if i == 1
            FCSizing = 0;
            dFC = 100;
        else
            FCSizing = 0;
            dFC = FCSizingMatrix(j,i-1) + 100;
        end
        [LCOEtemp,~,~,~,~,~,~,~,~,~,~] = RenewablesWithFuelCell(PercSolar(j),Scale(i),A,FCSizing);
        check = 0;
        while check == 0
             FCSizing = FCSizing + dFC;
            [LCOEtemp2,Penetrationtemp,~,~,~,~,~,~,~,~,~] = RenewablesWithFuelCell(PercSolar(j),Scale(i),A,FCSizing);
            if LCOEtemp2 > LCOEtemp
                FCSizing = FCSizing - dFC;
                dFC = dFC/5;
                if abs(LCOEtemp2 - LCOEtemp)/LCOEtemp < tol
                    FCSizing = FCSizing + 3*dFC;
                    check = 1;
                    disp(i);
                end
%             elseif Penetrationtemp == 1-A
%                 FCSizing = FCSizing - dFC;
%                 dFC = dFC/3;
%                 if dFC < 1
%                     check = 1;
%                     disp(i);
%                 end
            end
            FCSizing2 = FCSizing*EffCurve.MaxFCPower;
            FCSizingMatrix(j,i) = FCSizing2;
            [LCOEtemp,PenetrationTemp,WindTemp,SolarTemp,RoundTripEfficiencyTemp,StorageTemp,MasHydTemp,CurtRateTemp,NoStoragePen,MHTM,LCOEHTM] = RenewablesWithFuelCell(PercSolar(j),Scale(i),A,FCSizing);
            LCOEMatrix(j,i) = LCOEtemp;
            PenetrationMatrix(j,i) = PenetrationTemp;
            iWindMatrix(j,i) = WindTemp;
            iSolarMatrix(j,i) = SolarTemp;
            MasHydMatrix(j,i) = MasHydTemp;
%             if FCSizing < 1
%                 RoundTripEfficiencyTemp = 0;
%             end
            RTEMatrix(j,i) = RoundTripEfficiencyTemp;
            StorageMatrix(j,i) = StorageTemp;
            CurtRateMatrix(j,i) = CurtRateTemp;
            NoStoragePenMatrix(j,i) = NoStoragePen;
            MHTMMatrix(j,i) = MHTM;
            LCOEHTMMatrix(j,i) = LCOEHTM;
        end
    end
    toc
end


PenetrationWind = PenetrationMatrix(1,:);
FCSizingWind = FCSizingMatrix(1,:);
LCOEWind = LCOEMatrix(1,:);
iWindWind = iWindMatrix(1,:);
iSolarWind = iSolarMatrix(1,:);
MasHydWind = MasHydMatrix(1,:);
RTEWind = RTEMatrix(1,:);
StorageWind = StorageMatrix(1,:);
CurtRateWind = CurtRateMatrix(1,:);
NoStoragePenWind = NoStoragePenMatrix(1,:);
MHTMWind = MHTMMatrix(1,:);
LCOEHTMWind = LCOEHTMMatrix(1,:);

n = find(PenetrationWind == 1-A,1,'first');
LCOEWind((n+1):end) = [];
PenetrationWind((n+1):end) = [];
FCSizingWind((n+1):end) = [];
iWindWind((n+1):end) = [];
iSolarWind((n+1):end) = [];
MasHydWind((n+1):end) = [];
RTEWind((n+1):end) = [];
StorageWind((n+1):end) = [];
CurtRateWind((n+1):end) = [];
NoStoragePenWind((n+1):end) = [];
MHTMWind((n+1):end) = [];
LCOEHTMWind((n+1):end) = [];


for i = 1:imax
    if i == 1
        LCOEOptimalA(i) = min(LCOEMatrix(:,1));
        n = find(LCOEMatrix(:,1) == LCOEOptimalA(i));
        PenetrationOptimalA(i) = PenetrationMatrix(n,1);
        iWindOptimalA(i) = iWindMatrix(n,1);
        iSolarOptimalA(i) = iSolarMatrix(n,1);
        PercSolarOptimalA(i) = PercSolar(n);
        FCSizeOptimalA(i) = FCSizingMatrix(n,1);
        StorageOptimalA(i) = StorageMatrix(n,1);
        RTEOptimalA(i) = RTEMatrix(n,1);
        MasHydOptimalA(i) = MasHydMatrix(n,1);
        CurtRateOptimalA(i) = CurtRateMatrix(n,1);
        NSPOptimalA(i) = NoStoragePenMatrix(n,1);
        MHTMOptimal(i) = MHTMMatrix(n,1);
        LCOEHTMOptimal(i) = LCOEHTMMatrix(n,1);
    else
        dLCOEdPen = (LCOEMatrix(:,i) - LCOEOptimalA(i-1))./(PenetrationMatrix(:,i) - PenetrationOptimalA(i-1));
        dLCOEdPen(isnan(dLCOEdPen)) = 10*max(dLCOEdPen);
        dLCOEdPen(PenetrationMatrix(:,i) < PenetrationOptimalA(i-1)) = 10*max(dLCOEdPen);
        [M,I] = min(dLCOEdPen);
        LCOEOptimalA(i) = LCOEMatrix(I,i);
        PenetrationOptimalA(i) = PenetrationMatrix(I,i);
        iWindOptimalA(i) = iWindMatrix(I,i);
        iSolarOptimalA(i) = iSolarMatrix(I,i);
        PercSolarOptimalA(i) = PercSolar(I);
        FCSizeOptimalA(i) = FCSizingMatrix(I,i);
        StorageOptimalA(i) = StorageMatrix(I,i);
        RTEOptimalA(i) = RTEMatrix(I,i);
        MasHydOptimalA(i) = MasHydMatrix(I,i);
        CurtRateOptimalA(i) = CurtRateMatrix(I,i);
        NSPOptimalA(i) = NoStoragePenMatrix(I,i);
        MHTMOptimal(i) = MHTMMatrix(I,i);
        LCOEHTMOptimal(i) = LCOEHTMMatrix(I,i);
    end
end

n = find(PenetrationOptimalA == 1-A,1,'first');
LCOEOptimalA((n+1):end) = [];
PenetrationOptimalA((n+1):end) = [];
iWindOptimalA((n+1):end) = [];
iSolarOptimalA((n+1):end) = [];
PercSolarOptimalA((n+1):end) = [];
FCSizeOptimalA((n+1):end) = [];
StorageOptimalA((n+1):end) = [];
RTEOptimalA((n+1):end) = [];
MasHydOptimalA((n+1):end) = [];
CurtRateOptimalA((n+1):end) = [];
NSPOptimalA((n+1):end) = [];
MHTMOptimal((n+1):end) = [];
LCOEHTMOptimal((n+1):end) = [];



for i = 1:imax
    tic
    Scale(i) = maxScale*(i/imax);
    if i == 1
        FCSizing = 0;
        dFC = 100;
    else
        FCSizing = 0;
        dFC = FCSizingMatrix(j,i-1) + 100;
    end
    [LCOEtemp,~,~,~,~,~,~,~,~,~,~] = RenewablesWithFuelCell(1,Scale(i),A,FCSizing);
    check = 0;
    while check == 0
        FCSizing = FCSizing + dFC;
        [LCOEtemp2,Penetrationtemp,~,~,~,~,~,~,~,~,~] = RenewablesWithFuelCell(1,Scale(i),A,FCSizing);
        if LCOEtemp2 > LCOEtemp
            FCSizing = FCSizing - dFC;
            dFC = dFC/3;
            if abs(LCOEtemp2 - LCOEtemp)/LCOEtemp < tol
                FCSizing = FCSizing + 3*dFC;
                check = 1;
            end
%         elseif Penetrationtemp == 1-A
%             FCSizing = FCSizing - dFC;
%             dFC = dFC/3;
%             if dFC < 1
%                 check = 1;
%             end
        end
    end
    FCSizing2 = FCSizing*EffCurve.MaxFCPower;
    FCSizingSolar(i) = FCSizing2;
    [LCOEtemp,PenetrationTemp,WindTemp,SolarTemp,RoundTripEfficiencyTemp,StorageTemp,MasHydTemp,CurtRateTemp,NoStoragePen,MHTM,LCOEHTM] = RenewablesWithFuelCell(1,Scale(i),A,FCSizing);
    LCOESolar(i) = LCOEtemp;
    PenetrationSolar(i) = PenetrationTemp;
    iWindSolar(i) = WindTemp;
    iSolarSolar(i) = SolarTemp;
%     if FCSizing < 1
%         RoundTripEfficiencyTemp = 0;
%     end
    RTESolar(i) = RoundTripEfficiencyTemp;
    StorageSolar(i) = StorageTemp;
    MasHydSolar(i) = MasHydTemp;
    CurtRateSolar(i) = CurtRateTemp;
    NoStoragePenSolar(i) = NoStoragePen;
    MHTMSolar(i) = MHTM;
    LCOEHTMSolar(i) = LCOEHTM;
    toc
end

n = find(PenetrationSolar == 1-A,1,'first');
LCOESolar((n+1):end) = [];
PenetrationSolar((n+1):end) = [];
iWindSolar((n+1):end) = [];
iSolarSolar((n+1):end) = [];
MasHydSolar((n+1):end) = [];
RTESolar((n+1):end) = [];
StorageSolar((n+1):end) = [];
CurtRateSolar((n+1):end) = [];
NoStoragePenSolar((n+1):end) = [];
FCSizingSolar((n+1):end) = [];
MHTMSolar((n+1):end) = [];
LCOEHTMSolar((n+1):end) = [];




PercSolarMix = 0.5;
for i = 1:imax
    tic
    Scale(i) = maxScale*(i/imax);
    if i == 1
        FCSizing = 0;
        dFC = 100;
    else
        FCSizing = 0;
        dFC = FCSizingMatrix(j,i-1) + 100;
    end
    [LCOEtemp,~,~,~,~,~,~,~,~,~,~] = RenewablesWithFuelCell(PercSolarMix,Scale(i),A,FCSizing);
    check = 0;
    while check == 0
        FCSizing = FCSizing + dFC;
        [LCOEtemp2,Penetrationtemp,~,~,~,~,~,~,~,~,~] = RenewablesWithFuelCell(PercSolarMix,Scale(i),A,FCSizing);
        if LCOEtemp2 > LCOEtemp
            FCSizing = FCSizing - dFC;
            dFC = dFC/3;
            if abs(LCOEtemp2 - LCOEtemp)/LCOEtemp < tol
                FCSizing = FCSizing + 3*dFC;
                check = 1;
            end
%         elseif Penetrationtemp == 1-A
%             FCSizing = FCSizing - dFC;
%             dFC = dFC/3;
%             if dFC < 1
%                 check = 1;
%             end
        end
    end
    FCSizing2 = FCSizing*EffCurve.MaxFCPower;
    FCSizingMix(i) = FCSizing2;
    [LCOEtemp,PenetrationTemp,WindTemp,SolarTemp,RoundTripEfficiencyTemp,StorageTemp,MasHydTemp,CurtRateTemp,NoStoragePen,MHTM,LCOEHTM] = RenewablesWithFuelCell(PercSolarMix,Scale(i),A,FCSizing);
    LCOEMix(i) = LCOEtemp;
    PenetrationMix(i) = PenetrationTemp;
    iWindMix(i) = WindTemp;
    iSolarMix(i) = SolarTemp;
%     if FCSizing < 1
%         RoundTripEfficiencyTemp = 0;
%     end
    RTEMix(i) = RoundTripEfficiencyTemp;
    StorageMix(i) = StorageTemp;
    MasHydMix(i) = MasHydTemp;
    CurtRateMix(i) = CurtRateTemp;
    NoStoragePenMix(i) = NoStoragePen;
    MHTMMix(i) = MHTM;
    LCOEHTMMix(i) = LCOEHTM;
    toc
end

n = find(PenetrationMix == 1-A,1,'first');
LCOEMix((n+1):end) = [];
PenetrationMix((n+1):end) = [];
iWindMix((n+1):end) = [];
iSolarMix((n+1):end) = [];
MasHydMix((n+1):end) = [];
RTEMix((n+1):end) = [];
StorageMix((n+1):end) = [];
CurtRateMix((n+1):end) = [];
NoStoragePenMix((n+1):end) = [];
FCSizingMix((n+1):end) = [];
MHTMMix((n+1):end) = [];
LCOEHTMMix((n+1):end) = [];


reSOFCData.LCOE = LCOEOptimalA;
reSOFCData.Pen = PenetrationOptimalA;
reSOFCData.PercSolar = PercSolarOptimalA;
save('reSOFCData.mat','reSOFCData');

save('reSOFCRawData.mat');


figure(1)
plot(PenetrationWind,LCOEWind,'-b','LineWidth',2)
hold on
plot(PenetrationSolar,LCOESolar,'-r','LineWidth',2);
plot(PenetrationMix,LCOEMix,'-m','LineWidth',2);
plot(PenetrationOptimalA,LCOEOptimalA,'--xk','LineWidth',2);
hold off
xlabel('Annual Renewable Energy Penetration');
ylabel('Net LCOE ($/MWh)');
legend('100% Wind',' 100% Solar','50% Solar, 50% Wind','Optimal Mix');

figure(2)
plot(PenetrationOptimalA,PercSolarOptimalA,'--xk','LineWidth',2);
xlabel('Annual Renewable Energy Penetration');
ylabel('Installed Solar Capacity Fraction');

figure(3)
plot(PenetrationOptimalA,(iWindOptimalA/(0.32*365*24)),'--xb','LineWidth',2);
hold on
plot(PenetrationOptimalA,(iSolarOptimalA/(0.18*365*24)),'--xr','LineWidth',2);
hold off
%title('Optimal Installation by Type');
xlabel('Annual Renewable Energy Penetration');
ylabel('Installed Capacity (MW)');
legend('Wind','Solar');

FCPower = EffCurve.MaxFCPower;

figure(4)
plot(PenetrationWind,FCPower*FCSizingWind,'-b','LineWidth',2);
hold on
plot(PenetrationSolar,FCPower*FCSizingSolar,'-r','LineWidth',2);
plot(PenetrationMix,FCPower*FCSizingMix,'-m','LineWidth',2);
plot(PenetrationOptimalA,FCPower*FCSizeOptimalA,'--xk','LineWidth',2);
hold off
xlabel('Annual Renewable Energy Penetration');
ylabel('Installed Capacity (MW)');
legend('Wind','Solar','50% Wind, 50% Solar','Optimal');

figure(5)
plot(PenetrationWind,StorageWind,'-b','LineWidth',2);
hold on
plot(PenetrationSolar,StorageSolar,'-r','LineWidth',2);
plot(PenetrationMix,StorageMix,'-m','LineWidth',2);
plot(PenetrationOptimalA,StorageOptimalA,'--xk','LineWidth',2);
hold off
%title('Required Methane Storage Capcity');
xlabel('Annual Renewable Energy Penetration');
ylabel('Mols CH4');
legend('Wind','Solar','50% Wind, 50% Solar','Optimal');

n = find(RTEWind >0,1,'first');
RTEWind(1:(n-1)) = [];
PW = PenetrationWind(n:end);
if isempty(PW)
    PW = 0*RTEWind;
end


n = find(RTESolar >0,1,'first');
RTESolar(1:(n-1)) = [];
PS = PenetrationSolar(n:end);
if isempty(PS)
    PS = 0*RTESolar;
end

n = find(RTEMix >0,1,'first');
RTEMix(1:(n-1)) = [];
PM = PenetrationMix(n:end);
if isempty(PM)
    PM = 0*RTEMix;
end


n = find(RTEOptimalA >0,1,'first');
RTEOptimalA(1:(n-1)) = [];
PO = PenetrationOptimalA(n:end);
if isempty(PO)
    PO = 0*RTEOptimalA;
end

figure(6)
plot(PW,RTEWind,'-b','LineWidth',2);
hold on
plot(PS,RTESolar,'-r','LineWidth',2);
plot(PM,RTEMix,'-m','LineWidth',2);
plot(PO,RTEOptimalA,'--xk','LineWidth',2);
hold off
xlabel('Annual Renewable Energy Penetration');
ylabel('Roundtrip Electrical Efficiency');
legend('Wind','Solar','50% Wind, 50% Solar','Optimal');

%mas hyd from electrolysis
figure(7)
plot(PenetrationWind,MasHydWind,'-b','LineWidth',2);
hold on
plot(PenetrationSolar,MasHydSolar,'-r','LineWidth',2);
plot(PenetrationMix,MasHydMix,'-m','LineWidth',2);
plot(PenetrationOptimalA,MasHydOptimalA,'--xk','LineWidth',2);
hold off
xlabel('Annual Renewable Energy Penetration');
ylabel('H2 Production (Kg/Year)');
legend('Wind','Solar','50% Wind, 50% Solar','Optimal');

%mas hyd from HTM
figure(8)
plot(PenetrationWind,MHTMWind,'-b','LineWidth',2);
hold on
plot(PenetrationSolar,MHTMSolar,'-r','LineWidth',2);
plot(PenetrationMix,MHTMMix,'-m','LineWidth',2);
plot(PenetrationOptimalA,MHTMOptimal,'--xk','LineWidth',2);
hold off
xlabel('Annual Renewable Energy Penetration');
ylabel('H2 Production (Kg/Year)');
legend('Wind','Solar','50% Wind, 50% Solar','Optimal');

%yearly revenue from Hydrogen
figure(9)
plot(PenetrationWind,LCOEHTMWind,'-b','LineWidth',2);
hold on
plot(PenetrationSolar,LCOEHTMSolar,'-r','LineWidth',2);
plot(PenetrationMix,LCOEHTMMix,'-m','LineWidth',2);
plot(PenetrationOptimalA,LCOEHTMOptimal,'--xk','LineWidth',2);
hold off
xlabel('Annual Renewable Energy Penetration');
ylabel('$/Year');
legend('Wind','Solar','50% Wind, 50% Solar','Optimal');




%curtailment rate
figure(10)
plot(PenetrationWind,CurtRateWind,'-b','LineWidth',2);
hold on
plot(PenetrationSolar,CurtRateSolar,'-r','LineWidth',2);
plot(PenetrationMix,CurtRateMix,'-m','LineWidth',2);
plot(PenetrationOptimalA,CurtRateOptimalA,'--xk','LineWidth',2);
hold off
xlabel('Annual Renewable Energy Penetration');
ylabel('Curtailment Rate');
legend('Wind','Solar','50% Wind, 50% Solar','Optimal');

%penetration from renewables
figure(11)
plot(PenetrationWind,NoStoragePenWind,'-b','LineWidth',2);
hold on
plot(PenetrationSolar,NoStoragePenSolar,'-r','LineWidth',2);
plot(PenetrationMix,NoStoragePenMix,'-m','LineWidth',2);
plot(PenetrationOptimalA,NSPOptimalA,'--xk','LineWidth',2);
hold off
xlabel('Annual Renewable Energy Penetration');
ylabel('Penetration Without Energy Storage');
legend('Wind','Solar','50% Wind, 50% Solar','Optimal');

%penetration from storage
figure(12)
plot(PenetrationWind,PenetrationWind-NoStoragePenWind,'-b','LineWidth',2);
hold on
plot(PenetrationSolar,PenetrationSolar-NoStoragePenSolar,'-r','LineWidth',2);
plot(PenetrationMix,PenetrationMix-NoStoragePenMix,'-m','LineWidth',2);
plot(PenetrationOptimalA,PenetrationOptimalA-NSPOptimalA,'--xk','LineWidth',2);
hold off
xlabel('Annual Renewable Energy Penetration');
ylabel('Penetration From Energy Storage');
legend('Wind','Solar','50% Wind, 50% Solar','Optimal');


reSOFCData.LCOE = LCOEOptimalA;
reSOFCData.Pen = PenetrationOptimalA;
reSOFCData.PercSolar = PercSolarOptimalA;
save('reSOFCData.mat','reSOFCData');
end