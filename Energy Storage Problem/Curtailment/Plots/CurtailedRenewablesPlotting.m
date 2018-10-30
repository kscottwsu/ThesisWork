function [ ] = CurtailedRenewablesPlotting()
A = 0.5;%percentage of demand met by baseline suppy. max value = 0.707
imax = 200;
jmax = 200;
maxScale = 2;


for i = 1:imax
    tic
    Scale(i) = maxScale*(i/imax);
    [LCOEWind(i),PenetrationWind(i),~,~,CurtailmentRateWind(i)] = CurtailedRenewables(0,Scale(i),A);   
    toc
end


for i = 1:imax
    tic
    Scale(i) = maxScale*(i/imax);
    [LCOESolar(i),PenetrationSolar(i),~,~,CurtailmentRateSolar(i)] = CurtailedRenewables(1,Scale(i),A);   
    toc
end


for i = 1:imax
    tic
    Scale(i) = maxScale*(i/imax);
    [LCOEMix(i),PenetrationMix(i),~,~,CurtailmentRateMix(i)] = CurtailedRenewables(0.5,Scale(i),A);   
    toc
end


for j = 1:jmax
    tic
    PercSolar(j) = 0.4*((j-1)/(jmax-1));
    for i = 1:imax
        Scale(i) = maxScale*(i/imax);
        [LCOEtemp,PenetrationTemp,WindTemp,SolarTemp,CurtailmentRateTemp] = CurtailedRenewables(PercSolar(j),Scale(i),A);   
        LCOEMatrix(j,i) = LCOEtemp;
        PenetrationMatrix(j,i) = PenetrationTemp;
        iWindMatrix(j,i) = WindTemp;
        iSolarMatrix(j,i) = SolarTemp;
        CurtailmentMatrix(j,i) = CurtailmentRateTemp;
    end
    toc
end


for i = 1:imax
    if i == 1
        LCOEOptimalA(i) = min(LCOEMatrix(:,1));
        n = find(LCOEMatrix(:,1) == LCOEOptimalA(i));
        PenetrationOptimalA(i) = PenetrationMatrix(n,1);
        iWindOptimalA(i) = iWindMatrix(n,1);
        iSolarOptimalA(i) = iSolarMatrix(n,1);
        PercSolarOptimalA(i) = PercSolar(n);
        CurtailmentOptimalA(i) = CurtailmentMatrix(n,1);
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
        CurtailmentOptimalA(i) = CurtailmentMatrix(I,i);
    end
end

n = find(PenetrationOptimalA == 1-A,1,'first');
LCOEOptimalA((n+1):end) = [];
PenetrationOptimalA((n+1):end) = [];
iWindOptimalA((n+1):end) = [];
iSolarOptimalA((n+1):end) = [];
PercSolarOptimalA((n+1):end) = [];
CurtailmentOptimalA((n+1):end) = [];

save('rawCurtailmentData.mat');


figure(1)
plot(PenetrationWind,LCOEWind,'LineWidth',2);
hold on
plot(PenetrationSolar,LCOESolar,'LineWidth',2);
plot(PenetrationMix,LCOEMix,'LineWidth',2);
plot(PenetrationOptimalA,LCOEOptimalA,'LineWidth',2);
hold off
%title('Net LCOE for Different Mixes of Wind and Solar Energy');
xlabel('Annual Renewable Energy Penetration');
ylabel('Net LCOE ($/MWh)');
legend('100% Wind',' 100% Solar','50% Wind, 50% Solar','Optimal Mix');
ylim([0 250]);



figure(2)
plot(PenetrationWind,CurtailmentRateWind,'LineWidth',2);
hold on
plot(PenetrationSolar,CurtailmentRateSolar,'LineWidth',2);
plot(PenetrationMix,CurtailmentRateMix,'LineWidth',2);
plot(PenetrationOptimalA,CurtailmentOptimalA,'LineWidth',2);
hold off
%title('Curtailment Rates for Different Mixes of Wind and Solar Energy');
xlabel('Annual Renewable Energy Penetration');
ylabel('Curtailment Rate');
legend('100% Wind',' 100% Solar','50% Wind, 50% Solar','Optimal Mix');

figure(3)
plot(PenetrationOptimalA,PercSolarOptimalA,'LineWidth',2);
%title('Percentage of Solar Energy at Optimal Cost');
xlabel('Annual Renewable Energy Penetration');
ylabel('Installed Solar Capacity Fraction');

figure(4)
plot(PenetrationOptimalA,(iWindOptimalA/(0.32*365*24)),'LineWidth',2);
hold on
plot(PenetrationOptimalA,(iSolarOptimalA/(0.18*365*24)),'LineWidth',2);
hold off
%title('Installed Capacity of Solar and Wind at Optimal Cost')
xlabel('Annual Renewable Energy Penetration');
ylabel('Installed Capacity (MWh)');
legend('Installed Wind','Installed Solar');

CurtData.LCOE = LCOEOptimalA;
CurtData.Pen = PenetrationOptimalA;
CurtData.PercSolar = PercSolarOptimalA;
save('CurtData.mat','CurtData');
end