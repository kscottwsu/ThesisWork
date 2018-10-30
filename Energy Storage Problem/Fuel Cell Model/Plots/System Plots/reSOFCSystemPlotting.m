function [ ] = reSOFCSystemPlotting( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
imax = 200;
bDisp = 0;
bHyd = 0;
PrEC = 2000;

FlowMin = 5e-7;
FlowMax = 8e-6;
SteamUtilization = 0.77;
for i = 1:imax
    tic
    Flow(i) = FlowMin + (FlowMax-FlowMin)*(i-1)/(imax-1);    
    [EIn,FlowOut,SysEff,j,Vfc,CellEff,PElectroTemp,OxiEffTemp] = ElectrolysisCellSystemModel(Flow(i),SteamUtilization,PrEC,bHyd,bDisp);
    EC.SysEff(i) = SysEff;
    EC.current(i) = -j;
    EC.Voltage(i) = Vfc;
    EC.Flow(i) = FlowOut.CH4;
    EC.SysPower(i) = -1e3*EIn;
    EC.CellEff(i) = CellEff;
    EC.CellPower(i) = -1e3*PElectroTemp;
    EC.NetCellEff(i) = OxiEffTemp;
    disp(i/imax);
    toc
end

bHyd = 1;
for i = 1:imax
    tic
    Flow(i) = FlowMin + (FlowMax-FlowMin)*(i-1)/(imax-1);    
    [EIn,FlowOut,SysEff,j,Vfc,CellEff,PElectroTemp,OxiEffTemp] = ElectrolysisCellSystemModel(Flow(i),SteamUtilization,PrEC,bHyd,bDisp);
    ECH.SysEff(i) = SysEff;
    ECH.current(i) = -j;
    ECH.Voltage(i) = Vfc;
    ECH.Hyd(i) = FlowOut.H2;
    ECH.SysPower(i) = -1e3*EIn;
    ECH.CellEff(i) = CellEff;
    ECH.CellPower(i) = -1e3*PElectroTemp;
    ECH.NetCellEff(i) = OxiEffTemp;
    disp(i/imax);
    toc
end

bHyd = 0;
PrFC = 2000;
DS2CRatio = 3;
PHTM = 7;%Above 120 results in no hydrogen
FlowMax = 2e-6;
FlowMin = 1e-7;
FCUtilization = 0.7;
for i = 1:imax
    tic
    Flow(i) = FlowMin + (FlowMax - FlowMin)*(i-1)/(imax-1);
    [PowerOut,SysEff,HydProduction,j,Vfc,FCEff,PFCT,NetSysEffTemp,OxiEffTemp] = FuelCellSystemModel(Flow(i),DS2CRatio, PHTM,PrFC,bHyd,FCUtilization,bDisp);
    FC.SysEff(i) = SysEff;
    FC.current(i) = j;
    FC.Voltage(i) = Vfc;
    FC.SysPower(i) = 1e3*PowerOut;
    FC.Hyd(i) = HydProduction;
    FC.CellEff(i) = FCEff;
    FC.CellPower(i) = 1e3*PFCT;
    FC.NetSysEff(i) = NetSysEffTemp;
    FC.NetCellEff(i) = OxiEffTemp;
    disp(i/imax);
    toc
end

bHyd = 1;
for i = 1:imax
    tic
    Flow(i) = FlowMin + (FlowMax - FlowMin)*(i-1)/(imax-1);
    [PowerOut,SysEff,HydProduction,j,Vfc,FCEff,PFCT,NetSysEffTemp,OxiEffTemp] = FuelCellSystemModel(Flow(i),DS2CRatio, PHTM,PrFC,bHyd,FCUtilization,bDisp);
    FCH.SysEff(i) = SysEff;
    FCH.current(i) = j;
    FCH.Voltage(i) = Vfc;
    FCH.SysPower(i) = 1e3*PowerOut;
    FCH.Hyd(i) = HydProduction;
    FCH.CellEff(i) = FCEff;
    FCH.CellPower(i) = 1e3*PFCT;
    FCH.NetSysEff(i) = NetSysEffTemp;
    FCH.NetCellEff(i) = OxiEffTemp;
    disp(i/imax);
    toc
end


temp = 0:0.001:1;

%% Cell Properties in System Model
%V vs I plus Cell Power
figure(1)
yyaxis left
plot(FC.current,FC.Voltage,'-b','LineWidth',2);
hold on
plot(EC.current,EC.Voltage,'-b','LineWidth',2);
plot(0*temp,2*max(EC.Voltage)*temp,'-k','LineWidth',2);
hold off
xlabel('Current Density (A/cm^2)');
ylabel('Voltage (V)');
xlim([(min(EC.current)-0.1) (max(FC.current)+0.1)])
ylim([(min(FC.Voltage)-0.1) (max(EC.Voltage)+0.1)])
yyaxis right
plot(FC.current,FC.CellPower,'-r','LineWidth',2);
hold on
plot(EC.current,EC.CellPower,'-r','LineWidth',2);
plot(0*temp,2*max(FC.CellPower)*temp,'-k','LineWidth',2)
hold off
ylabel('Power Density (W/cm^2)');
ylim([(min(EC.CellPower)-0.1) (max(FC.CellPower)+0.1)])

%% System Net Efficiency

figure(3)
% yyaxis left
plot(EC.SysPower,EC.SysEff,'-b','LineWidth',2);
hold on
plot(FC.SysPower,FC.SysEff,'-b','LineWidth',2);
plot(FC.SysPower,FC.NetSysEff,'b--','LineWidth',2);%includes hydrogen production by HTM
plot(0*temp,2*max(FC.NetSysEff)*temp,'-k','LineWidth',2);
hold off
ylabel('System Efficiency');
ylim([(min(FC.SysEff)-0.1) (max(FC.NetSysEff)+0.1)])
% yyaxis right
% plot(FC.SysPower,FC.Hyd,'-r','LineWidth',2);
% hold on
% plot(0*temp,2*max(FC.Hyd)*temp,'-k','LineWidth',2);
% hold off
% ylabel('H2 Production (Mols/sec cm^2)');
xlabel('System Power Density (W/cm^2)');
xlim([(min(EC.SysPower)-0.1) (max(FC.SysPower)+0.1)])
% ylim([(0) (1.1*max(FC.Hyd))]);

%% H2 Production From Electrolysis Mode



figure(9)
% yyaxis left
plot(EC.SysPower,EC.SysEff,'-b','LineWidth',2);
hold on
plot(FC.SysPower,FC.SysEff,'-b','LineWidth',2);
plot(0*temp,2*max(FC.SysEff)*temp,'-k','LineWidth',2);
plot(FCH.SysPower,FCH.SysEff,'-r','LineWidth',2);
plot(ECH.SysPower,ECH.SysEff,'-r','LineWidth',2);
hold off
ylabel('System Efficiency');
xlabel('System Power Density (W/cm^2)');
ylim([0.4 1.05]);
legend('Methane-based System','Methane-based System','Hydrogen-based System','Hydrogen-based System');
%need to add legend

figure(10)
plot(EC.SysPower,EC.CellEff,'-b','LineWidth',2);
hold on
plot(FC.SysPower,FC.CellEff,'-b','LineWidth',2);
plot(0*temp,2*max(FC.CellEff)*temp,'-k','LineWidth',2);
plot(FCH.SysPower,FCH.CellEff,'-r','LineWidth',2);
plot(ECH.SysPower,ECH.CellEff,'-r','LineWidth',2);
hold off
ylabel('Cell Efficiency');
xlabel('System Power Density (W/cm^2)');
ylim([0.4 1.05]);
legend('Methane-based System','Methane-based System','Hydrogen-based System','Hydrogen-based System');
%need to add legend



for i = 1:length(EC.current)
    EC.SysEff2(i) = EC.SysEff(end+1-i);
    EC.SysPower2(i) = EC.SysPower(end+1-i);
end

reSOFC.Eff = [EC.SysEff2,EC.SysEff2(end),FC.SysEff(1),FC.SysEff];
reSOFC.HydProductionVec = [0*EC.SysEff2,0,1e-6*FC.Hyd(1),FC.Hyd];
reSOFC.SysPower = [EC.SysPower2,1e-6*EC.SysPower2(end),1e-6*FC.SysPower(1),FC.SysPower];

reSOFC.ECHSysPower = [0,ECH.SysPower];
reSOFC.HydECEff = [ECH.SysEff(1),ECH.SysEff];
reSOFC.HydEC = [0, ECH.Hyd];

EffCurve.Eff =  reSOFC.Eff;
EffCurve.HydProduction = reSOFC.HydProductionVec;
EffCurve.Power = reSOFC.SysPower;
EffCurve.MaxFCPower = 0.6;
EffCurve.MinFCPower = -1.85;
EffCurve.ElectrolysisModeEff = reSOFC.HydECEff;
EffCurve.ElectrolysisModePower = reSOFC.ECHSysPower;
EffCurve.ElectrolysisModeHyd = reSOFC.HydEC;

save('FC.mat','FC');
save('EC.mat','EC');
save('ECH.mat','ECH');
save('EffCurve.mat','EffCurve');
end