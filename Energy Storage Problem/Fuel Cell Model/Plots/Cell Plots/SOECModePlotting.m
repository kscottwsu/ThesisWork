function [ ] = SOECModePlotting(bHyd,PrEC,f)
imax = 200;
bDisp = 0;

Flow = 3e-6;
SUMax = 0.995;
SUMin = 0.01;

if bHyd == 0
    SUMax = 0.77;
end

for i = 1:imax
    tic
    SteamUtilization(i) = SUMin + (SUMax-SUMin)*(i-1)/(imax-1);
    [EIn,FlowOut,SysEff,j,Vfc,ECEff,PElectroTemp,OxiEffTemp ] = ElectrolysisCellSystemModel(Flow,SteamUtilization(i),PrEC,bHyd,bDisp);
    EffES(i) = SysEff;
    currentE(i) = j;
    VoltageE(i) = Vfc;
    FlowE(i) = FlowOut.CH4;
    PECS(i) = 1e3*EIn;
    EffEC(i) = ECEff;
    PEC(i) = 1e3*PElectroTemp;
    NEffE(i) = OxiEffTemp;
    disp(i/imax);
    toc
end

figure(13)
hold on
if bHyd == 1
    plot(currentE,VoltageE,'LineWidth',2);
else
    plot(currentE,VoltageE,'m','LineWidth',2);
end

figure(14)
hold on
if bHyd == 1
    plot(currentE,PEC,'LineWidth',2);
else
    plot(currentE,PEC,'m','LineWidth',2);
end

figure(15)
hold on
if bHyd == 1
    plot(currentE,NEffE,'LineWidth',2);
elseif bHyd ~= 1
    plot(currentE,EffEC,'m','LineWidth',2);
    plot(currentE,NEffE,'--m','LineWidth',2);
end
% 
% bDisp = 1;
% SteamUtilization = 0.77;
% Flow = 3e-6;
% [EIn,FlowOut,SysEff,j,Vfc,SysEff,PElectroTemp,OxiEffTemp ] = ElectrolysisCellSystemModel(Flow,SteamUtilization,PrEC,bHyd,bDisp);
end

