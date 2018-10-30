function   LCOEComparison( )
load('BatteryData.mat');
load('CurtData.mat');
load('reSOFCData.mat');


figure(1)
hold off
plot(CurtData.Pen,CurtData.LCOE,'-','LineWidth',2);
hold on
plot(BatteryData.Pen,BatteryData.LCOE,'-','LineWidth',2);
plot(reSOFCData.Pen,reSOFCData.LCOE,'-','LineWidth',2);
xlabel('Annual Renewable Energy Penetration');
ylabel('Net LCOE ($/MWh)');
legend('Curtailment','Battery Energy Storage','reSOFC-based Energy Storage');


figure(2)
hold off
plot(CurtData.Pen,CurtData.PercSolar,'--x','LineWidth',2);
hold on
plot(BatteryData.Pen,BatteryData.PercSolar,'--x','LineWidth',2);
plot(reSOFCData.Pen,reSOFCData.PercSolar,'--x','LineWidth',2);
legend('Curtailment','Battery Energy Storage','reSOFC-based Energy Storage');
xlabel('Annual Renewable Energy Penetration');
ylabel('Installed Solar Capacity Fraction');
end

