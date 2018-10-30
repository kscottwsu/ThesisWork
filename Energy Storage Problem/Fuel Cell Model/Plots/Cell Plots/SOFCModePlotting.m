function [ ] = SOFCModePlotting(bHyd,PrFC,f)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
imax =  200;
DS2CRatio = 3;
PHTM = 5;

bDisp = 0;

Flow = 1e-6;
FCUMax = 0.995;
FCUMin = 0.01;

for i = 1:imax
    tic
    FCUtilization(i) = FCUMin + (FCUMax-FCUMin)*(i-1)/(imax-1);
    [PowerOut,SysEff,HydProduction,j,Vfc,FCEff,PFCT,NetSysEffTemp,OxiEffTemp] = FuelCellSystemModel(Flow,DS2CRatio, PHTM,PrFC,bHyd,FCUtilization(i),bDisp);
    EffFCS(i) = SysEff;
    currentFC(i) = j;
    VoltageFC(i) = Vfc;
    if i > 1
        if VoltageFC(i) > VoltageFC(i-1)
            disp('error');
        end
    end
         
    PFCS(i) = 1e3*PowerOut;
    HydFC(i) = HydProduction;
    EffFC(i) = FCEff;
    PFC(i) = 1e3*PFCT;
    NEffFCS(i) = NetSysEffTemp;
    NCEffFC(i) = OxiEffTemp;
    disp(i/imax);
    toc
end

figure(10)
hold on
if bHyd == 1
    plot(currentFC,VoltageFC,'LineWidth',2);
else
    plot(currentFC,VoltageFC,'m','LineWidth',2);
end

figure(11)
hold on
if bHyd == 1
    plot(currentFC,PFC,'LineWidth',2);
else
    plot(currentFC,PFC,'m','LineWidth',2);
end

figure(12)
hold on
if bHyd == 1
    plot(currentFC,EffFC,'LineWidth',2);
else
    plot(currentFC,EffFC,'m','LineWidth',2);
    plot(currentFC,NCEffFC,'m--','LineWidth',2);
end

% bDisp = 1;
% FCUtilization = 0.95;
% Flow = 1e-6;
% [PowerOut,SysEff,HydProduction,j,Vfc,FCEff,PFCT,NetSysEffTemp,OxiEffTemp] = FuelCellSystemModel(Flow,DS2CRatio, PHTM,PrFC,bHyd,FCUtilization,bDisp);
end

