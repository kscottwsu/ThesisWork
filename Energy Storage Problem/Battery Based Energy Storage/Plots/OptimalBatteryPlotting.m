function [ ] = OptimalBatteryPlotting( )
%Calculates Optimal LCOE with lithium ion battery based energy storage
A = 0.5;%percentage of demand met by baseline suppy. max value = 0.707
imax = 150;
jmax = 150;
maxScale = 3.6;

for j = 1:jmax
    tic
    PercSolar(j) = ((j-1)/(jmax-1));
    for i = 1:imax
        Scale(i) = maxScale*(i/imax);
        if i == 1
            BatterySizing = 0;
            dBattery = 100;
        else
            BatterySizing = 0;
            dBattery = BatterySizingMatrix(j,i-1) + 100;
        end

        [LCOEtemp,Penetrationtemp,~,~] = RenewablesWithBattery(PercSolar(j),Scale(i),A,BatterySizing);
        check = 0;
        tol = 1e-6;
        while check == 0
            BatterySizing = BatterySizing + dBattery;
            [LCOEtemp2,Penetrationtemp2,~,~] = RenewablesWithBattery(PercSolar(j),Scale(i),A,BatterySizing);
            if (LCOEtemp2-LCOEtemp)/(Penetrationtemp2-Penetrationtemp) > 0
                BatterySizing = BatterySizing - dBattery; 
                dBattery = dBattery/5;
                if abs(LCOEtemp2-LCOEtemp)/LCOEtemp < tol
                    BatterySizing = BatterySizing + 5*dBattery;
                    check = 1;
                    disp(i);
                end
            end
        end
        BatterySizingMatrix(j,i) = BatterySizing;
        [LCOEtemp,PenetrationTemp,WindTemp,SolarTemp] = RenewablesWithBattery(PercSolar(j),Scale(i),A,BatterySizing);
        LCOEMatrix(j,i) = LCOEtemp;
        PenetrationMatrix(j,i) = PenetrationTemp;
        iWindMatrix(j,i) = WindTemp;
        iSolarMatrix(j,i) = SolarTemp;
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
        BatterySizeOptimalA(i) = BatterySizingMatrix(n,1);
    else
        dLCOEdPen = (LCOEMatrix(:,i) - LCOEOptimalA(i-1))./(PenetrationMatrix(:,i) - PenetrationOptimalA(i-1));
        %dLCOEdPen = dLCOEdPen + 1e-12;
        dLCOEdPen(dLCOEdPen <0) = 10*max(dLCOEdPen);
        dLCOEdPen(isnan(dLCOEdPen)) = 10*max(dLCOEdPen);
        dLCOEdPen(PenetrationMatrix(:,i) < PenetrationOptimalA(i-1)) = 10*max(dLCOEdPen);
        [M,I] = min(dLCOEdPen);
        LCOEOptimalA(i) = LCOEMatrix(I,i);
        PenetrationOptimalA(i) = PenetrationMatrix(I,i);
        iWindOptimalA(i) = iWindMatrix(I,i);
        iSolarOptimalA(i) = iSolarMatrix(I,i);
        PercSolarOptimalA(i) = PercSolar(I);
        BatterySizeOptimalA(i) = BatterySizingMatrix(I,i);
    end
end


n = find(PenetrationOptimalA >= 1-A,1,'first');
LCOEOptimalA((n+1):end) = [];
PenetrationOptimalA((n+1):end) = [];
iWindOptimalA((n+1):end) = [];
iSolarOptimalA((n+1):end) = [];
PercSolarOptimalA((n+1):end) = [];
BatterySizeOptimalA((n+1):end) = [];

save('rawBatteryData.mat');





figure(1)
plot(PenetrationMatrix(1,:),LCOEMatrix(1,:),'LineWidth',2)
hold on
plot(PenetrationMatrix(jmax,:),LCOEMatrix(jmax,:),'LineWidth',2);
plot(PenetrationMatrix(round(jmax/2),:),LCOEMatrix(round(jmax/2),:),'LineWidth',2);
plot(PenetrationOptimalA,LCOEOptimalA,'-x','LineWidth',2);
hold off
%title('Net LCOE For Mixes of Wind and Solar Energy');
xlabel('Annual Renewable Energy Penetration');
ylabel('Net LCOE ($/MWh)');
legend('100% Wind',' 100% Solar','50% Solar, 50% Wind','Optimal Mix');

figure(2)
plot(PenetrationOptimalA,(iWindOptimalA/(0.32*365*24)),'-x','LineWidth',2);
hold on
plot(PenetrationOptimalA,(iSolarOptimalA/(0.18*365*24)),'-x','LineWidth',2);
hold off
%title('Optimal Installation by Type');
xlabel('Annual Renewable Energy Penetration');
ylabel('Installed Capacity (MW)');
legend('Wind','Solar');

figure(3)
plot(PenetrationOptimalA,PercSolarOptimalA,'-x','LineWidth',2);
%title('Optimal Percentage of Energy Coming From Solar');
xlabel('Annual Renewable Energy Penetration');
ylabel('Installed Solar Capacity Fraction');


figure(4)
plot(PenetrationMatrix(1,:),BatterySizingMatrix(1,:),'LineWidth',2);
hold on
plot(PenetrationMatrix(jmax,:),BatterySizingMatrix(jmax,:),'LineWidth',2);
plot(PenetrationMatrix(round(jmax/2),:),BatterySizingMatrix(round(jmax/2),:),'LineWidth',2);
plot(PenetrationOptimalA,BatterySizeOptimalA,'-x','LineWidth',2);
hold off
%title('Optimal Battery Sizing');
xlabel('Annual Renewable Energy Penetration');
ylabel('Installed Capacity (MWh)');
legend('Wind','Solar','50% Wind, 50% Solar','Optimal');

BatteryData.LCOE = LCOEOptimalA;
BatteryData.Pen = PenetrationOptimalA;
BatteryData.PercSolar = PercSolarOptimalA;
save('BatteryData.mat','BatteryData');
end