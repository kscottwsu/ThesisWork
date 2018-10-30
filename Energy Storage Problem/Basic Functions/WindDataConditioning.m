function WindDataConditioning( ~ )
global StorageEffSOFC StorageEffSOEC Methanation HydrogenEff CurtailmentPerc



%BatteryAnalysis(0.2);




% 
% Eefficiency2 = 0.35;
% StorageEffSOEC = 0.90;
% StorageEffSOFC = Eefficiency2/StorageEffSOEC;
% Methanation = 0.98;
% HydrogenEff = 0.02;
% CurtailmentPerc = 30;
% 
% CurtailmentPerc = 0;
% [data,~] = GetData(0.18);
% 
% for i = 1:100
%     x(i) = i;
%     temp = -data(data(:,2)<0,2);
%     temp2 = data(data(:,2)>0,2);
%     y(i) = prctile(temp,i);
%     y2(i) = prctile(temp2,i);
% end
% 
% figure(99)
% hold off
% plot(x,y);
% hold on
% plot(x,y2);
% ylabel('Normalized Power');
% xlabel('Demand Percentile');

% 
% Eefficiency2 = 0.55;
% StorageEffSOEC = 0.90;
% StorageEffSOFC = Eefficiency2/StorageEffSOEC;
% Methanation = 0.98;
% HydrogenEff = 0.02;
% CurtailmentPerc = 30;
% 
% CurtailmentPerc = 0;
% [data,~] = GetData(0.18);
% 
% for i = 1:100
%     x(i) = i;
%     temp = -data(data(:,2)<0,2);
%     temp2 = data(data(:,2)>0,2);
%     y(i) = prctile(temp,i);
%     y2(i) = prctile(temp2,i);
% end
% 
% figure(99)
% hold on
% plot(x,y);
% hold on
% plot(x,y2);
% 
% Eefficiency2 = 0.75;
% StorageEffSOEC = 0.90;
% StorageEffSOFC = Eefficiency2/StorageEffSOEC;
% Methanation = 0.98;
% HydrogenEff = 0.02;
% CurtailmentPerc = 30;
% 
% CurtailmentPerc = 0;
% [data,~] = GetData(0.18);
% 
% for i = 1:100
%     x(i) = i;
%     temp = -data(data(:,2)<0,2);
%     temp2 = data(data(:,2)>0,2);
%     y(i) = prctile(temp,i);
%     y2(i) = prctile(temp2,i);
% end
% 
% figure(99)
% hold on
% plot(x,y);
% hold on
% plot(x,y2);
% legend('35% Efficiency: Electrolyzer','35% Efficiency: Fuel Cell','55% Efficiency: Electrolyzer','55% Efficiency: Fuel Cell','75% Efficiency: Electrolyzer','75% Efficiency: Fuel Cell');



% 
% 
% Eefficiency2 = 0.55;
% StorageEffSOEC = 0.90;
% StorageEffSOFC = Eefficiency2/StorageEffSOEC;
% Methanation = 0.98;
% HydrogenEff = 0.02;
% 
% CurtailmentPerc = 0;
% [data,analysis1] = GetData(0.18);
% Demand1 = data(:,2)*analysis1.RFCsize;
% 
% 
% 
% CurtailmentPerc = 20;
% [data,analysis2] = GetData(0.18);
% Demand2 = data(:,2)*analysis2.RFCsize;
% time = data(:,1);
% 
% figure(1)
% hold off
% plot(time,Demand1);
% hold on
% plot(time,Demand2);
% xlabel('Time (s)');
% ylabel('Demand');
% legend('0% psuedo-curtailment','20% psuedo-curtailment');











%     
% MethanationFraction = zeros(1000,1);
% for i = 1:1000
%     SOECeff(i) = i/1000;
%     [MethanationFraction(i)] = ThermalBalance(SOECeff(i));
% end
% 
% MethanationFraction(MethanationFraction > 1) = 1;
% 
% 
% plot(SOECeff,MethanationFraction);
% ylabel('Thermodynamically Stable Methanation Fraction');
% xlabel('SOEC Electrical Efficiency');
% axis([0 1 0 1+0.1])
    
    Eefficiency2 = 0.55;
    StorageEffSOEC = 0.90;
    StorageEffSOFC = Eefficiency2/StorageEffSOEC;
    Methanation = 0.98;
    HydrogenEff = 0.10;
    CurtailmentPerc = 10;


    PercSolar = zeros(1,201);
    UtilizationBoth = PercSolar;
    StorageFractionBoth = PercSolar;
    RFCcost = PercSolar;
    SolarCost = PercSolar;
    WindCost = PercSolar;
    HydrogenDollars = PercSolar;
    CurtailmentLossPerc = PercSolar;
    StorageCost = PercSolar;
    RFCutil = PercSolar;
    Storage = zeros(105404,200);
    Demand = Storage;
    MaxRampRate = PercSolar;
    RampRate = Storage;
    
for i = 0:200
    [data,analysis] = GetData(i/200);
    PercSolar(i+1) = analysis.PercSolar;
    UtilizationBoth(i+1) = analysis.UtilizationBoth;
    StorageFractionBoth(i+1) = analysis.StorageFractionBoth;
    %TotalCost(i+1) = analysis{i+1}.TotalCost;
    StorageCost(i+1) = analysis.StorageCost;
    RFCcost(i+1) = analysis.RFCcost;
    RFCutil(i+1) = analysis.FCUtil;
    SolarCost(i+1) = analysis.SolarCost;
    WindCost(i+1) = analysis.WindCost;
    Storage(:,i+1) = data(:,3);
    Demand(:,i+1) = data(:,2);
    HydrogenDollars(i+1) = analysis.EffCostReduction;
    CurtailmentLossPerc(i+1) = analysis.CurtailmentLossDemand;
    time = data(:,1);
    RampRate(:,i+1) = data(:,4);
    MaxRampRate(i+1) = analysis.MaxRampRate;
end
CostperMWh = SolarCost + WindCost + RFCcost + StorageCost;

minUtil = min(UtilizationBoth);
nMinUtil = find(UtilizationBoth == minUtil);

minCost = min(CostperMWh-HydrogenDollars);
nMinCost = find(CostperMWh-HydrogenDollars == minCost);


figure(1)
yyaxis left
plot(PercSolar,UtilizationBoth);
ylabel('Utilization Fraction');
yyaxis right
plot(PercSolar,StorageFractionBoth);
ylabel('Max Storage Fraction');

figure(2)
plot(time,Storage(:,1));
hold on
plot(time,Storage(:,nMinUtil));
plot(time,Storage(:,nMinCost));
plot(time,Storage(:,end));
ylabel('Storage Fraction');
xlabel('Time (s)');
legend('Wind','Minimum Utilization','Minimum Cost','Solar');
hold off

figure(3)
plot(time,RampRate(:,1));
hold on
plot(time,RampRate(:,nMinUtil));
plot(time,RampRate(:,nMinCost));
plot(time,RampRate(:,end));
ylabel('Ramp Rate %/Min')
legend('Wind','Minimum Utilization','Minimum Cost','Solar');

figure(7)
plot(PercSolar,RFCcost);
hold on
plot(PercSolar,WindCost);
plot(PercSolar,SolarCost);
plot(PercSolar,HydrogenDollars);
plot(PercSolar,StorageCost);
hold off
ylabel('$/MWh');
xlabel('Solar Fraction');
legend('Reversible Fuel Cell','Wind','Solar','Hydrogen Offset','Storage Costs');

figure(8)
plot(PercSolar,CostperMWh-HydrogenDollars);
ylabel('Effective $/MWh');
xlabel('PercSolar');

figure(9)
yyaxis left
plot(PercSolar,CurtailmentLossPerc);
xlabel('Solar Fraction');
ylabel('Curtailment Loss %');
yyaxis right
plot(PercSolar,MaxRampRate)
ylabel('Maximum Ramp Rate %/Min ');


figure(10)
plot(PercSolar,RFCutil);
xlabel('Solar Fraction');
ylabel('RFC Capacity Factor');

data = []; %#ok<NASGU>
analysis = [];%#ok<NASGU>
PercSolar = [];%#ok<NASGU>
UtilizationBoth = [];%#ok<NASGU>
StorageFractionBoth = [];%#ok<NASGU>
CostperMWh = [];%#ok<NASGU>
RFCcost = [];%#ok<NASGU>
SolarCost = [];%#ok<NASGU>
WindCost = [];%#ok<NASGU>
Storage = [];%#ok<NASGU>
Demand = [];%#ok<NASGU>
HydrogenDollars = [];%#ok<NASGU>
time = [];%#ok<NASGU>

PercCurt = zeros(1,200);
PercSolar = PercCurt;
UtilizationBoth = PercSolar;
StorageFractionBoth = PercSolar;
RFCcost = PercSolar;
StorageCost = PercSolar;
SolarCost = PercSolar;
RFCutil = PercSolar;
WindCost = PercSolar;
HydrogenDollars = PercSolar;
CurtailmentLossPerc = PercSolar;
Storage = zeros(105404,200);
Demand = Storage;
MaxRampRate = PercSolar;
RampRate = Storage;



Eefficiency2 = 0.35;
StorageEffSOEC = 0.95;
StorageEffSOFC = Eefficiency2/StorageEffSOEC;
Methanation = 0.98;
HydrogenEff = 0.02;
CurtailmentPerc = 0;

for i = 1:200
    PercCurt(i) = 90*(i-1)/200;
    CurtailmentPerc = PercCurt(i);
    [data,analysis] = GetData(0.2);
    %PercSolar(i+1) = analysis.PercSolar;
    UtilizationBoth(i) = analysis.UtilizationBoth;
    StorageFractionBoth(i) = analysis.StorageFractionBoth;
    %TotalCost(i+1) = analysis{i+1}.TotalCost;
    StorageCost(i) = analysis.StorageCost;
    RFCcost(i) = analysis.RFCcost;
    RFCutil(i) = analysis.FCUtil;
    SolarCost(i) = analysis.SolarCost;
    WindCost(i) = analysis.WindCost;
    Storage(:,i) = data(:,3);
    Demand(:,i) = data(:,2);
    HydrogenDollars(i) = analysis.EffCostReduction;
    CurtailmentLossPerc(i) = analysis.CurtailmentLossDemand;
    time = data(:,1);
    RampRate(:,i) = data(:,4);
    %Remainder(:,i) = data(:,5);
    MaxRampRate(i) = analysis.MaxRampRate;
end


CostperMWh = SolarCost + WindCost + RFCcost + StorageCost;

minUtil = min(UtilizationBoth);
nMinUtil = find(UtilizationBoth == minUtil);

minCost = min(CostperMWh-HydrogenDollars);
nMinCost = find(CostperMWh-HydrogenDollars == minCost);


figure(1)
yyaxis left
plot(PercCurt,UtilizationBoth);
ylabel('Utilization Fraction');
yyaxis right
plot(PercCurt,StorageFractionBoth);
ylabel('Max Storage Fraction');


figure(7)
plot(PercCurt,RFCcost);
hold on
plot(PercCurt,WindCost);
plot(PercCurt,SolarCost);
plot(PercCurt,HydrogenDollars);
plot(PercCurt,StorageCost);
hold off
ylabel('$/MWh');
xlabel('Solar Fraction');
legend('Reversible Fuel Cell','Wind','Solar','Hydrogen Offset','Storage Costs');

figure(8)
plot(PercCurt,CostperMWh-HydrogenDollars);
ylabel('Effective $/MWh');
xlabel('Curtailment Percentile');

figure(9)
yyaxis left
plot(PercCurt,CurtailmentLossPerc);
xlabel('Curtailment Percentile');
ylabel('Curtailment Loss %');
yyaxis right
plot(PercCurt,MaxRampRate)
ylabel('Maximum Ramp Rate %/Min ');



figure(10)
plot(PercCurt,RFCutil);
xlabel('Pseudo-Curtailment Percentage');
ylabel('RFC Capacity Factor');




data = []; %#ok<NASGU>
analysis = [];%#ok<NASGU>
PercSolar = [];%#ok<NASGU>
UtilizationBoth = [];%#ok<NASGU>
StorageFractionBoth = [];%#ok<NASGU>
CostperMWh = [];%#ok<NASGU>
RFCcost = [];%#ok<NASGU>
SolarCost = [];%#ok<NASGU>
WindCost = [];%#ok<NASGU>
Storage = [];%#ok<NASGU>
Demand = [];%#ok<NASGU>
HydrogenDollars = [];%#ok<NASGU>
time = [];%#ok<NASGU>

PercCurt = zeros(1,200);
%Remainder = PercCurt;
PercSolar = PercCurt;
UtilizationBoth = PercSolar;
StorageFractionBoth = PercSolar;
RFCcost = PercSolar;
StorageCost = PercSolar;
RFCutil = PercCurt;
SolarCost = PercSolar;
WindCost = PercSolar;
HydrogenDollars = PercSolar;
CurtailmentFraction = PercSolar;
Storage = zeros(105404,200);
Demand = Storage;
MaxRampRate = PercSolar;
RampRate = Storage;



Eefficiency2 = 0.55;
StorageEffSOEC = 0.95;
StorageEffSOFC = Eefficiency2/StorageEffSOEC;
Methanation = 0.98;
HydrogenEff = 0.02;
CurtailmentPerc = 2;

for i = 1:200
    PercCurt(i) = 90*(i-1)/200;
    CurtailmentPerc = PercCurt(i);
    [data,analysis] = GetData(0.18);
    %PercSolar(i+1) = analysis.PercSolar;
    UtilizationBoth(i) = analysis.UtilizationBoth;
    StorageFractionBoth(i) = analysis.StorageFractionBoth;
    %TotalCost(i+1) = analysis{i+1}.TotalCost;
    RFCcost(i) = analysis.RFCcost;
    StorageCost(i) = analysis.StorageCost;
    RFCutil(i) = analysis.FCUtil;
    SolarCost(i) = analysis.SolarCost;
    WindCost(i) = analysis.WindCost;
    Storage(:,i) = data(:,3);
    Demand(:,i) = data(:,2);
    %Remainder(:,i) = data(:,5);
    HydrogenDollars(i) = analysis.EffCostReduction;
    CurtailmentLossPerc(i) = analysis.CurtailmentLossDemand;
    time = data(:,1);
    RampRate(:,i) = data(:,4);
    MaxRampRate(i) = analysis.MaxRampRate;
end


CostperMWh = SolarCost + WindCost + RFCcost + StorageCost;

minUtil = min(UtilizationBoth);
nMinUtil = find(UtilizationBoth == minUtil);

minCost = min(CostperMWh-HydrogenDollars);
nMinCost = find(CostperMWh-HydrogenDollars == minCost);


figure(1)
yyaxis left
plot(PercCurt,UtilizationBoth);
ylabel('Utilization Fraction');
yyaxis right
plot(PercCurt,StorageFractionBoth);
ylabel('Max Storage Fraction');


figure(7)
plot(PercCurt,RFCcost);
hold on
plot(PercCurt,WindCost);
plot(PercCurt,SolarCost);
plot(PercCurt,HydrogenDollars);
plot(PercCurt,StorageCost);
hold off
ylabel('$/MWh');
xlabel('Solar Fraction');
legend('Reversible Fuel Cell','Wind','Solar','Hydrogen Offset','Storage Costs');

figure(8)
plot(PercCurt,CostperMWh-HydrogenDollars);
ylabel('Effective $/MWh');
xlabel('Curtailment Percentile');

figure(9)
yyaxis left
plot(PercCurt,CurtailmentLossPerc);
xlabel('Solar Fraction');
ylabel('Curtailment Loss %');
yyaxis right
plot(PercCurt,MaxRampRate)
ylabel('Maximum Ramp Rate %/Min ');







Eefficiency2 = 0.75;
StorageEffSOEC = 0.95;
StorageEffSOFC = Eefficiency2/StorageEffSOEC;
Methanation = 0.98;
HydrogenEff = 0.02;
CurtailmentPerc = 2;

for i = 1:200
    PercCurt(i) = 90*(i-1)/200;
    CurtailmentPerc = PercCurt(i);
    [data,analysis] = GetData(0.18);
    %PercSolar(i+1) = analysis.PercSolar;
    UtilizationBoth(i) = analysis.UtilizationBoth;
    StorageFractionBoth(i) = analysis.StorageFractionBoth;
    %TotalCost(i+1) = analysis{i+1}.TotalCost;
    RFCcost(i) = analysis.RFCcost;
    StorageCost(i) = analysis.StorageCost;
    RFCutil(i) = analysis.FCUtil;
    SolarCost(i) = analysis.SolarCost;
    WindCost(i) = analysis.WindCost;
    Storage(:,i) = data(:,3);
    Demand(:,i) = data(:,2);
    HydrogenDollars(i) = analysis.EffCostReduction;
    CurtailmentLossPerc(i) = analysis.CurtailmentLossDemand;
    time = data(:,1);
    %Remainder(:,i) = data(:,5);
    RampRate(:,i) = data(:,4);
    MaxRampRate(i) = analysis.MaxRampRate;
end


CostperMWh = SolarCost + WindCost + RFCcost + StorageCost;

minUtil = min(UtilizationBoth);
nMinUtil = find(UtilizationBoth == minUtil);

minCost = min(CostperMWh-HydrogenDollars);
nMinCost = find(CostperMWh-HydrogenDollars == minCost);


figure(1)
yyaxis left
plot(PercCurt,UtilizationBoth);
ylabel('Utilization Fraction');
yyaxis right
plot(PercCurt,StorageFractionBoth);
ylabel('Max Storage Fraction');


figure(7)
plot(PercCurt,RFCcost);
hold on
plot(PercCurt,WindCost);
plot(PercCurt,SolarCost);
plot(PercCurt,HydrogenDollars);
plot(PercCurt,StorageCost);
hold off
ylabel('$/MWh');
xlabel('Solar Fraction');
legend('Reversible Fuel Cell','Wind','Solar','Hydrogen Offset','Storage Costs');

figure(8)
plot(PercCurt,CostperMWh-HydrogenDollars);
ylabel('Effective $/MWh');
xlabel('Curtailment Percentile');

figure(9)
yyaxis left
plot(PercCurt,CurtailmentLossPerc);
xlabel('Solar Fraction');
ylabel('Curtailment Loss %');
yyaxis right
plot(PercCurt,MaxRampRate)
ylabel('Maximum Ramp Rate %/Min ');









Eefficiency2 = 1;
StorageEffSOEC = 1;
StorageEffSOFC = Eefficiency2/StorageEffSOEC;
Methanation = 0.98;
HydrogenEff = 0;
CurtailmentPerc = 0;

for i = 1:200
    PercCurt(i) = 90*(i-1)/200;
    CurtailmentPerc = PercCurt(i);
    [data,analysis] = GetData(0.18);
    %PercSolar(i+1) = analysis.PercSolar;
    UtilizationBoth(i) = analysis.UtilizationBoth;
    StorageFractionBoth(i) = analysis.StorageFractionBoth;
    %TotalCost(i+1) = analysis{i+1}.TotalCost;
    RFCcost(i) = analysis.RFCcost;
    StorageCost(i) = analysis.StorageCost;
    RFCutil(i) = analysis.FCUtil;
    SolarCost(i) = analysis.SolarCost;
    WindCost(i) = analysis.WindCost;
    Storage(:,i) = data(:,3);
    Demand(:,i) = data(:,2);
    HydrogenDollars(i) = analysis.EffCostReduction;
    CurtailmentLossPerc(i) = analysis.CurtailmentLossDemand;
    time = data(:,1);
    %Remainder(:,i) = data(:,5);
    RampRate(:,i) = data(:,4);
    MaxRampRate(i) = analysis.MaxRampRate;
end


CostperMWh = SolarCost + WindCost + RFCcost + StorageCost;

minUtil = min(UtilizationBoth);
nMinUtil = find(UtilizationBoth == minUtil);

minCost = min(CostperMWh-HydrogenDollars);
nMinCost = find(CostperMWh-HydrogenDollars == minCost);


figure(1)
yyaxis left
plot(PercCurt,UtilizationBoth);
ylabel('Utilization Fraction');
yyaxis right
plot(PercCurt,StorageFractionBoth);
ylabel('Max Storage Fraction');


figure(7)
plot(PercCurt,RFCcost);
hold on
plot(PercCurt,WindCost);
plot(PercCurt,SolarCost);
plot(PercCurt,HydrogenDollars);
plot(PercCurt,StorageCost);
hold off
ylabel('$/MWh');
xlabel('Solar Fraction');
legend('Reversible Fuel Cell','Wind','Solar','Hydrogen Offset','Storage Costs');

figure(8)
plot(PercCurt,CostperMWh-HydrogenDollars);
ylabel('Effective $/MWh');
xlabel('Curtailment Percentile');

figure(9)
yyaxis left
plot(PercCurt,CurtailmentLossPerc);
xlabel('Solar Fraction');
ylabel('Curtailment Loss %');
yyaxis right
plot(PercCurt,MaxRampRate)
ylabel('Maximum Ramp Rate %/Min ');

































Eefficiency2 = 0.75;
StorageEffSOEC = 0.95;
StorageEffSOFC = Eefficiency2/StorageEffSOEC;
Methanation = 0.98;
HydrogenEff = 0.02;
CurtailmentPerc = 30;

[data,~] = GetData(0.18);

for i = 1:100
    x(i) = i;
    temp = -data(data(:,2)<0,2);
    temp2 = data(data(:,2)>0,2);
    y(i) = prctile(temp,i);
    y2(i) = prctile(temp2,i);
end

figure(99)
hold on
plot(x,y);
hold on
plot(x,y2);






%[data,analysis] = GetData(0.3);















data = []; %#ok<NASGU>
analysis = [];%#ok<NASGU>
PercSolar = [];%#ok<NASGU>
UtilizationBoth = [];%#ok<NASGU>
StorageFractionBoth = [];%#ok<NASGU>
CostperMWh = [];%#ok<NASGU>
RFCcost = [];%#ok<NASGU>
SolarCost = [];%#ok<NASGU>
WindCost = [];%#ok<NASGU>
Storage = [];%#ok<NASGU>
Demand = [];%#ok<NASGU>
HydrogenDollars = [];%#ok<NASGU>
time = [];%#ok<NASGU>






kmax = 20;
jmax = 20;
imax = 20;

PercSolar = zeros(imax+1,1);
UtilizationBoth = PercSolar;
StorageFractionBoth = PercSolar;
RFCcost = PercSolar;
StorageCost = PercSolar;
SolarCost = PercSolar;
WindCost = PercSolar;
HydrogenDollars = PercSolar;
minCost = zeros(jmax,1);
nMinCost = minCost;
minStorage = minCost;
DHydrogen = minCost;
Eefficiency = minCost;
dMinCost = zeros(jmax-1,1);
dHdC = dMinCost;
CurFraction = minCost;
MaxRampRate = PercSolar;
MaxRRate = minCost;
CurtailmentP = zeros(kmax,1);
%find impact of optimal cost vs electrical efficiency
for k = 1:kmax
    CurtailmentP(k) = 40*(k-1)/kmax;
    CurtailmentPerc = CurtailmentP(k);
    for j = 1:jmax
        Eefficiency(j) = 0.3+0.65*(j-1)/jmax;
        StorageEffSOEC = 0.90;
        StorageEffSOFC = Eefficiency(j)/StorageEffSOEC;
        HydrogenEff = 0.65*1/jmax;
        
        for i = 0:imax
            PercSolar(i+1) = 0.08+0.16*i/imax;
            [~,analysis] = GetData(PercSolar(i+1));%scale it to solar fraction of 0.08 - 0.24
            %PercSolar(i+1) = analysis.PercSolar;
            UtilizationBoth(i+1) = analysis.UtilizationBoth;
            StorageFractionBoth(i+1) = analysis.StorageFractionBoth;
            StorageCost(i+1) = analysis.StorageCost;
            RFCcost(i+1) = analysis.RFCcost;
            SolarCost(i+1) = analysis.SolarCost;
            WindCost(i+1) = analysis.WindCost;
            HydrogenDollars(i+1) = analysis.EffCostReduction;
            CurtailmentLossPerc(i+1) = analysis.CurtailmentLossPerc;
            MaxRampRate(i+1) = analysis.MaxRampRate;
        end
        CostperMWh = SolarCost + WindCost + RFCcost + StorageCost;
        minCost(j,k) = min(CostperMWh);
        nMinCost(j,k) = find(CostperMWh == minCost(j,k));
        minStorage(j,k) = StorageFractionBoth(nMinCost(j,k));
        DHydrogen(j,k) = HydrogenDollars(nMinCost(j,k));%hydrogen $ generated at min cost(excluding hydrogen production)
        CurFraction(j,k) = CurtailmentFraction(nMinCost(j,k));
        MaxRRate(j,k) = MaxRampRate(nMinCost(j,k));
        if j > 1
            dMinCost(j,k) = minCost(j-1,k)-minCost(j,k);
            dHdC(j,k) = DHydrogen(j-1,k)/dMinCost(j,k);
        end
    end
end

figure(21)
yyaxis left
plot(Eefficiency,minCost);
ylabel('$/MWh');
xlabel('Electrical Efficiency');
yyaxis right
plot(Eefficiency,0.08+0.16*nMinCost/imax)
ylabel('Solar Fraction');

BreakEven = ones(length(dHdC),1);
figure(22)
plot(Eefficiency(1:end),dHdC);
hold on
plot(Eefficiency(1:end),BreakEven);
legend('Hydrogen Income/Cost','Break-even Point')
hold off

figure(23)
yyaxis left
plot(Eefficiency,minStorage);
ylabel('Storage Fraction');

%find impact of equivalent value of 1% hydrogen output vs electrical
%efficiency


figure(24)
yyaxis left
plot(Eefficiency,MaxRRate);
xlabel('Electrical Efficiency');
ylabel('Maximum Ramp Rate %/Min');
yyaxis right
plot(Eefficiency,CurFraction);
ylabel('Curtailment Loss Fraction');
xlabel('Electrical Efficiency');
end

function [Data,Analysis] = GetData(PercSolar)
global StorageEffSOFC StorageEffSOEC HydrogenEff CurtailmentPerc

StorageEff = StorageEffSOFC*StorageEffSOEC;
Analysis.PercSolar = PercSolar;
RawWindData = load('WindData.mat');
RawData = RawWindData.RawWindData.Spring;

n = 0:(length(RawData(:,1))-1);
t = n*5*60;%convert from 5 minutes to seconds
RawData(:,4) = t';%time in seconds

iDemand = 300/(60*60)*sum(RawData(:,2));%in MWh
iWind = 300/(60*60)*sum(RawData(:,1));%roughly 4782 MW installed capacity
iSolar = 300/(60*60)*sum(RawData(:,3));%solar data in MWh, roughly 25 MW peak power


CapacityFactor =  0.1531;%from nrel solar pv calc for Yakima Wa %0.1973;% for pheonix Az % 0.1451 for Pullman Wa
InstalledSolar = (iSolar/CapacityFactor)/(365*24); %go from MWh to MW installed



RawData(:,1) = (iDemand/iWind)*RawData(:,1);%in MW increase wind to match demand
RawData(:,3) = (iDemand/iSolar)*RawData(:,3);
RawData(:,9) = (1-PercSolar)*RawData(:,1) + PercSolar*RawData(:,3);
%c = iWind/(365*24*4782) = 0.2741;

InstalledWind = (1-PercSolar)*(iDemand/iWind)*4782;%(in MW)  4782 from BPA website
InstalledSolar = PercSolar*(iDemand/iSolar)*InstalledSolar;
iWind = iDemand*(1-PercSolar);
iSolar = iDemand*PercSolar;
if isnan(InstalledWind)
    InstalledWind = 0;
end



if isnan(InstalledSolar)
    InstalledSolar = 0;
end


RawData(:,5) = RawData(:,2) - RawData(:,1);%net power as a function of time
RawData(:,6) = RawData(:,2) - RawData(:,3);% in MW
RawData(:,7) = RawData(:,2) - RawData(:,9);

iStorageBoth = (300/(60*60))*sum((RawData(:,7) > 0).*RawData(:,7));%amount of storage needed to supply demand in MWh

iLossesBoth = (iStorageBoth/StorageEff - iStorageBoth);%Amount of additional energy needed by storage efficiency losses
iFractionBoth = iLossesBoth/iDemand;

RawData(:,9) = RawData(:,9)*(1+iFractionBoth);
RawData(:,7) = RawData(:,2) - RawData(:,9);%net power as a function of time

InstalledSolar = (1+iFractionBoth)*InstalledSolar;
InstalledWind = (1+iFractionBoth)*InstalledWind;
%costs from NREL Energy Technology Cost and Performance Data


iWind = iWind*(1+iFractionBoth);
iSolar = iSolar*(1+iFractionBoth);




%CostperMWhHydrogenTank = 14.69e3;%price per MWh, for hydrogen storage at 700 bar from DOE FY 2015 annual progress report 
for j = 1:3% converges towards proper energy balance
    iStorageBothT = 0*RawData(:,7);
    iHydrogen = iStorageBothT;
    iDemandSatisfied1 = 0;
    for i = 1:length(RawData(:,5))% energy storage in MWh
        if i > 1
            if RawData(i,7) > 0
                %iElectrolysis1 = iElectrolysis1;
                iDemandSatisfied1 = iDemandSatisfied1 + 300/(60*60)*RawData(i,7);
                iStorageBothT(i) = iStorageBothT(i-1) - 300/(60*60)*RawData(i,7)/StorageEffSOFC;
                iHydrogen(i) = iHydrogen(i-1);%iHydrogen(i-1) + 300/(60*60)*RawData(i,7)*HydrogenEff;
            else
                %iDemandSatisfied1 = iDemandSatisfied1;
                iElectrolysis1 = iElectrolysis1 - 300/(60*60)*RawData(i,7)*StorageEffSOEC;
                iStorageBothT(i) = iStorageBothT(i-1) - 300/(60*60)*RawData(i,7)*StorageEffSOEC;
                iHydrogen(i) = iHydrogen(i-1) -300/(60*60)*RawData(i,7)*HydrogenEff;
            end
        else
            if RawData(i,7) > 0
                iElectrolysis1 = 0;
                iDemandSatisfied1 = 300/(60*60)*RawData(i,7);
                iStorageBothT(i) =  -300/(60*60)*RawData(i,7)/StorageEffSOFC;
                iHydrogen(i) = 0;
            else
                iElectrolysis1 = -300/(60*60)*RawData(i,7)*StorageEffSOEC;
                iDemandSatisfied1 = 0;
                iStorageBothT(i) =  -300/(60*60)*RawData(i,7)*StorageEffSOEC;
                iHydrogen(i) = -300/(60*60)*RawData(i,7)*HydrogenEff;
            end
        end
    end

    iStorageBothT = (iStorageBothT - min(iStorageBothT))/(iDemand);
    EnergySurplus = (iStorageBothT(end) - iStorageBothT(1));
    RawData(:,9) = RawData(:,9)*(1-EnergySurplus*(2/3));
    RawData(:,7) = RawData(:,2) - RawData(:,9);
    InstalledSolar = InstalledSolar*(1-EnergySurplus*(2/3));
    InstalledWind = InstalledWind*(1-EnergySurplus*(2/3));
    iWind = iWind*(1-EnergySurplus*(2/3));
    iSolar = iSolar*(1-EnergySurplus*(2/3));
end



PmaxBoth = max(RawData(:,7));
PminBoth = -min(RawData(:,7));
%curtails to percentile to minimize reversible fuel cell sizing
temp =  100 - CurtailmentPerc;
if PmaxBoth > PminBoth
    CurtailmentThreshold = prctile(RawData((RawData(:,7)>0),7),temp);
    if CurtailmentThreshold > PminBoth
        PminBoth = CurtailmentThreshold;
        RawData((RawData(:,7)<PminBoth),7) = PminBoth;
    end
    RawData((RawData(:,7)>CurtailmentThreshold),7) = CurtailmentThreshold;
    Ptemp = CurtailmentThreshold;
else
    CurtailmentThreshold = -prctile(-RawData((RawData(:,7)<0),7),temp);
    if -CurtailmentThreshold < PmaxBoth
        PmaxBoth = -CurtailmentThreshold;
        RawData((RawData(:,7)>PmaxBoth),7) = PmaxBoth;
    end
    RawData((RawData(:,7)<CurtailmentThreshold),7) = CurtailmentThreshold;
    Ptemp = -CurtailmentThreshold;
end


for j = 1:5

    iStorageBothT = 0*RawData(:,7);
    iHydrogen = iStorageBothT;
    iDemandSatisfied2 = 0;
    for i = 1:length(RawData(:,5))% energy storage in MWh
        if i > 1
            if RawData(i,7) > 0
                %iElectrolysis2 = iElectrolysis2;
                iDemandSatisfied2 = iDemandSatisfied2 + 300/(60*60)*RawData(i,7); 
                iStorageBothT(i) = iStorageBothT(i-1) - 300/(60*60)*RawData(i,7)/StorageEffSOFC;
                iHydrogen(i) = iHydrogen(i-1);%iHydrogen(i-1) + 300/(60*60)*RawData(i,7)*HydrogenEff;
            else
                %iDemandSatisfied2 = iDemandSatisfied2;
                iElectrolysis2 = iElectrolysis2 - 300/(60*60)*RawData(i,7)*StorageEffSOEC;
                iStorageBothT(i) = iStorageBothT(i-1) - 300/(60*60)*RawData(i,7)*StorageEffSOEC;
                iHydrogen(i) = iHydrogen(i-1) -300/(60*60)*RawData(i,7)*HydrogenEff;
            end
        else
            if RawData(i,7) > 0
                iElectrolysis2 = 0;
                iDemandSatisfied2 = 300/(60*60)*RawData(i,7);
                iStorageBothT(i) =  -300/(60*60)*RawData(i,7)/StorageEffSOFC;
                iHydrogen(i) = 0;
            else
                iElectrolysis2 = -300/(60*60)*RawData(i,7)*StorageEffSOEC;
                iDemandSatisfied2 = 0;
                iStorageBothT(i) =  -300/(60*60)*RawData(i,7)*StorageEffSOEC;
                iHydrogen(i) = -300/(60*60)*RawData(i,7)*HydrogenEff;
            end
        end
    end
    iStorageBothT = (iStorageBothT - min(iStorageBothT))/(iDemand);

    CError = (iStorageBothT(1) - iStorageBothT(end))/(max(iStorageBothT));
    RawData(RawData(:,7)>0,7) = RawData(RawData(:,7)>0,7)*(1-CError/3);
end

Remainder = RawData(:,2) - RawData(:,9) - RawData(:,7);
SOECcurtailment = Remainder;%curtailment from SOEC sizing optimization
SOECcurtailment(SOECcurtailment>0) = 0;
SOECcurtailment = -SOECcurtailment;
Remainder(Remainder<0) = 0;
Ebudget = sum(Remainder);
guess = Ebudget/(length(Remainder));
Curtailment = 0*Remainder;%curtailment from constant supply(power not supplied by the energy storage system must be constant throughout the year) constraint

CError = 10;
while abs(CError) >= 1e-6
    Etotal = 0;
    FuelCellMode = RawData(:,7) + Remainder;
    FuelCellMode(FuelCellMode<0) = 0;
    
    for i = 1:length(FuelCellMode)
        if FuelCellMode(i) > guess
            Etotal = Etotal+guess;
            FuelCellMode(i) = FuelCellMode(i) - guess;
        else
            Curtailment(i) = (guess - FuelCellMode(i));
            Etotal = Etotal + FuelCellMode(i);
            FuelCellMode(i) = 0;
        end
    end
    CError = (Etotal - Ebudget)/Ebudget;
    guess = guess*(1-CError);
end

DemandRemainder = guess*ones(length(Curtailment),1);
Curtailment = Curtailment + SOECcurtailment;
TCurt = sum(Curtailment);
TRenewables = sum(RawData(:,9));

RenewablesCurtailmentPerc = 100*(1-(TRenewables-TCurt)/TRenewables);
RawData(RawData(:,7)>0,7) = FuelCellMode(RawData(:,7)>0);

%RawData(RawData(:,7)>Pbudget,7) = RawData(RawData(:,7)>Pbudget,7) - Pbudget*0;
%Remainder = RawData(:,2) - RawData(:,9) - RawData(:,7);
%Remainder(Remainder<0) = 0;
%RawData(:,7) = RawData(:,7)+Remainder;
%Remainder = RawData(:,2) - RawData(:,9) - RawData(:,7);
% Remainder(Remainder<0) = 0;
% Pconstant = mean(Remainder);
% EnonRenewable = Pconstant*RawData(end,4)/(60*60);%in MWh
% Remainder = Remainder - Pconstant;
% 
% demand = RawData(:,7);
% demand(demand<0) = 0;
% demand = demand+Remainder;%curtail SOFC if in SOFC mode
% demand(RawData(:,7)<0) = 0;%curtail renewables going directly to demand if in SOEC mode
% Imbalance = (1-StorageEff)*(300/(60*60))*sum(demand(demand<0));%in MWh
% RawData((RawData(:,7)>0),7) = demand((RawData(:,7)>0));
% 
% 
iStorageBothT = 0*RawData(:,7);
iHydrogen = iStorageBothT;
iDemandSatisfied2 = 0;
iElectrolysis2 = 0;
for i = 1:length(RawData(:,5))% energy storage in MWh
    if i > 1
        if RawData(i,7) > 0
            %iElectrolysis2 = iElectrolysis2;
            iDemandSatisfied2 = iDemandSatisfied2 + 300/(60*60)*RawData(i,7); 
            iStorageBothT(i) = iStorageBothT(i-1) - 300/(60*60)*RawData(i,7)/StorageEffSOFC;
            iHydrogen(i) = iHydrogen(i-1);%iHydrogen(i-1) + 300/(60*60)*RawData(i,7)*HydrogenEff;
        else
            %iDemandSatisfied2 = iDemandSatisfied2;
            iElectrolysis2 = iElectrolysis2 - 300/(60*60)*RawData(i,7)*StorageEffSOEC;
            iStorageBothT(i) = iStorageBothT(i-1) - 300/(60*60)*RawData(i,7)*StorageEffSOEC;
            iHydrogen(i) = iHydrogen(i-1) -300/(60*60)*RawData(i,7)*HydrogenEff;
        end
    else
        if RawData(i,7) > 0
            iElectrolysis2 = 0;
            iDemandSatisfied2 = 300/(60*60)*RawData(i,7);
            iStorageBothT(i) =  -300/(60*60)*RawData(i,7)/StorageEffSOFC;
            iHydrogen(i) = 0;
        else
            iElectrolysis2 = -300/(60*60)*RawData(i,7)*StorageEffSOEC;
            iDemandSatisfied2 = 0;
            iStorageBothT(i) =  -300/(60*60)*RawData(i,7)*StorageEffSOEC;
            iHydrogen(i) = -300/(60*60)*RawData(i,7)*HydrogenEff;
        end
    end
end
iStorageBothT = (iStorageBothT - min(iStorageBothT))/(iDemand);

CError2 = (iStorageBothT(1) - iStorageBothT(end))/(max(iStorageBothT));



%% ramp rate limitations/ battery storage integration
pmax = max(RawData(:,7));
pmin = -min(RawData(:,7));
maximumRampRate = 0.02*max(pmax,pmin);
maximumChange = 5*maximumRampRate;%time steps size is 5 minutes

Battery = 0*RawData(:,7);

for i = 2:length(RawData(:,7))
    if RawData(i,7) > RawData(i-1,7) + maximumChange;
        Battery(i) = RawData(i,7) - RawData(i-1,7) - maximumChange;
        RawData(i,7) = RawData(i-1,7) + maximumChange;
    elseif RawData(i,7) < RawData(i-1,7) - maximumChange;
        Battery(i) = RawData(i,7) -RawData(i-1,7) + maximumChange;
        RawData(i,7) = RawData(i-1,7) - maximumChange;
    else
        Battery(i) = 0;
    end
end







Pmin = -min(RawData(:,7));
Pmax = max(RawData(:,7));
Ptemp = max(Pmin,Pmax);
RawData(:,7) = RawData(:,7)/Ptemp;
Analysis.RFCsize = Ptemp;

Analysis.StorageFractionBoth = max(iStorageBothT);
Analysis.UtilizationBoth = iStorageBoth/iDemand;
iDemand2 = iDemand - (300/(60*60))*sum(DemandRemainder);
Analysis.CurtailmentLossDemand = 100*(iDemand-iDemand2)/iDemand;

%assumes that energy imbalance is dealt with through constant scaling
%instead of time dependent choice in when to scale output.  


RampRate = zeros(length(RawData(:,7)),1);
RampRate(2:end) = abs(100*(60)*(RawData(2:end,7) - RawData(1:end-1,7))./(RawData(2:end,4)-RawData(1:end-1,4)));%ramp rate in %/min
Analysis.MaxRampRate = max(RampRate);

%% storage cost

%CH4 storage
%10 million USD / 2*10e5 m^3 for salt cavern storage
%from: tightness and suitability evaluation of abandoned salt caverns ...
%260 KJ/mol
%1 MWh = 60*60 MJ
%assume pressure of 200 bar = 200*100 kPa
%assume temperature of 298 K
R = 8.314;%J*K-1 mol-1
T = 333.15;
p = 160*100*1000;%in pascals

Ech4 = Analysis.StorageFractionBoth*iDemand;%energy used to create methane in MWh
molCH4 = (Ech4*60*60)/0.260;
VolCH4 = 1.5*molCH4*R*T/p;%1/3 cushion gas
CostStorageCH4 = 10*10^6*VolCH4/(2*10^5);%upfront $ cost from CH4 storage 

%assume equal CO2 storage = maxCH4Storage - currentCH4 Storage, in order to
%preserve molar balance
molCO2 = molCH4;
VolCO2 = 1.5*molCO2*R*T/p;
CostStorageCO2 = 10*10^6*VolCO2/(2*10^5);


CostStorageO2 = CostStorageCO2;

molH2O = 4*molCH4;
molmassH2O = 18.01528/1000;%kg per mol
densityH2O = 998;%kg/m^3
VolumeH2O = (molH2O*molmassH2O)/densityH2O;%volume in m^3

CostStorageH2O = 50.4*VolumeH2O;

THydrogen = max(iHydrogen)*(60*60);%convert from MWh to MJ
THydrogen = ((THydrogen/0.23713)/2.01)/1000;%237.13 KJ/mol. 2.01 g/mol
CostHydrogen = THydrogen*2.00;%2.00,NREL production target according to Dustin,3.10 $/kg at delivery for NREL Target,  Yearly Income
%LCOE calculation assumes 30 year life time.



InstalledRFC = Ptemp;
CostperMWRSOFC = 1.5e6;% 1.5e3/kW 2020 goal from nrel Stationary Fuel Cell Evaluation in $/kW
RFCcost = InstalledRFC*CostperMWRSOFC;

n = 30;
DiscountRate = 0.056;% from eia cost determinations
DegredationRate = 0.01;%1% per year
OandM = 10;%assume roughly comparable to wind/solar power. $/MWh
TOandM = OandM*iDemand;


for i = 1:n
    if i == 1
        LCOE1 = RFCcost;
        LCOE2 = 0;
        LCOE2nodegrade = 0;
        LCOEStorage = CostStorageCH4 + CostStorageCO2 + CostStorageO2 +CostStorageH2O; 
        LCOEH = 0;
    end
    LCOE1 = LCOE1 + TOandM/(1+DiscountRate)^i;
    LCOE2 = LCOE2 + (iDemand2*(1-DegredationRate)^i)/((1+DiscountRate)^i);
    LCOE2nodegrade = LCOE2nodegrade + iDemand2/((1+DiscountRate)^i);
    LCOEH = LCOEH + (CostHydrogen*(1-DegredationRate)^i)/(1+DiscountRate)^i;%integrated demand including affects of degredation
end
LCOERFC = LCOE1/LCOE2;
HLCOE = LCOEH/LCOE2;
LCOEStorage = LCOEStorage/LCOE2nodegrade;

Analysis.RFCcost = LCOERFC;%this is a cost estimate because there is no lcoe data
Analysis.EffCostReduction = HLCOE;
Analysis.StorageCost = LCOEStorage;

CostperMWhSolarUnSub = 74.2;% from eia 2015 energy generation pdf for estimated lcoe for 2022
CostperMWhSolarSub = 58.2;%in $/MWh
altCostSolar = (CostperMWhSolarSub*iSolar/iDemand)*(iDemand/iDemand2);
Analysis.SolarCost = altCostSolar;
Analysis.SolarCostHigh = altCostSolar*0.26/0.195;%account for realistic capacity factor in Washington

CostperMWhWindUnSub = 58.5;
CostperMWhWindSub = 50.9;

altCostWind = (CostperMWhWindSub*iWind/iDemand)*(iDemand/iDemand2);

Analysis.WindCost = altCostWind;
Analysis.WindCostHigh = altCostWind*0.42/0.25;

FCUtilizationFactor = sum(abs(RawData(:,7)))/length(RawData(:,7));
Analysis.FCUtil = FCUtilizationFactor;


Data(:,1) = RawData(:,4);%time in seconds
Data(:,2) = RawData(:,7);%Net demand for entirely wind + solar System
Data(:,3) = iStorageBothT;
Data(:,4) = RampRate;
end

function [MethanationFraction] = ThermalBalance(EffSOEC)

%at 650 degrees C
%fuel cell mode
%H2 + 1/2 O2 -> H20, delh = -247 KJ/mol
%CH4 + H20 -> 3H2 + CO, delh = 224 KJ/mol
%H2O + CO -> H2 + CO2, delh = 36 KJ/mol

%CH4 +2 H20 -> 4 H2 + CO2, delh = 260 KJ/mol

%net balaced reaction FC mode
%CH4 + 2 O2 -> 2 H2O + CO2, delh = -988 EpowerKJ/mol, +260TpowerKJ/mol, net -728 KJ/mol
%Tpower used for thermal management
%always a balanced reaction ideally
% 
% FuelUtilization = 0.96;%CH4 utilization rate
% %remaining hydrogen fraction will be 4*remaining CH4
% 
% 
% EffSOFC = 0.6;
% Pout = 1;
% 
% EFC = 0.988/EffSOFC;%MJ/mol
% mols = (1/EFC)/FuelUtilization;%mols of CH4 consumed
% molsH = 4*mols*(1-FuelUtilization);%mols produced of H2
% Hydrogen = (molsH/2.01)/1000;%in kg
% 

%finds the thermodynamically neutral Reformation percentage that
%corresponds to a particular SOEC efficiency


%net balanced reaction electrolysis mode
%delh = 988 EpowerKJ/mol, - 260TpowerKJ/mol, net 728 KJ/mol
%Tpower used to increase efficiency


T = 650+273.15;%temp in kelvin
A = 30;%surface area m^2
ConvectiveHeatTransfer = 20e-6*(T-298)*A;%MW/M2 K
StefanBoltz = 5.670373e-14;%MW/(M2 K-4)
Emissivity = 0.1;%reflective metal
RadiativeHeatTransfer = StefanBoltz*Emissivity*A*(T^4-298^4);
Q = -ConvectiveHeatTransfer + -RadiativeHeatTransfer;

Hhydrogen = (0.988/EffSOEC)*(1-EffSOEC);
Hstorage = Hhydrogen + 0.260;




% 
% EStorage = (0.988/EffSOEC) - 0.260;%MJ/mol of CH4 stored
% HStorage = (1-EffSOEC)*EStorage;
% EHydrogen = (0.988/EffSOEC);%MJ/mol of 4 H2 Produced
% HHydrogen = (1-EffSOEC)*EHydrogen;


MethanationFraction = (-Q-Hhydrogen)/(Hstorage - Hhydrogen);
if MethanationFraction > 1
    MethanationFraction = 1;
elseif MethanationFraction < 0
    MethanationFraction = 0;
end
% 
% EffectiveEfficiency = (1*0.988)/(MethanationFraction*Hstorage
% 




end

function BatteryAnalysis (PercSolar)
%shows infeasiblity of seasonal shifting with batteries

Ebat = 0.95;%efficieincy in each direction
SelfDischargePM = 0.06;%percent per month
SelfDischarge = SelfDischargePM/(30*24*60*60);%percent per second, scales linearly with SOC
DegredationRate = 0.01;%percent per year
CostB = 100e3;%cost per MWh

StorageEff = Ebat^2;
Analysis.PercSolar = PercSolar;
RawWindData = load('WindData.mat');
RawData = RawWindData.RawWindData.Spring;

n = 0:(length(RawData(:,1))-1);
t = n*5*60;%convert from 5 minutes to seconds
RawData(:,4) = t';%time in seconds

iDemand = 300/(60*60)*sum(RawData(:,2));%in MWh
iWind = 300/(60*60)*sum(RawData(:,1));%roughly 4782 MW installed capacity
iSolar = 300/(60*60)*sum(RawData(:,3));%solar data in MWh, roughly 25 MW peak power


CapacityFactor =  0.1531;%from nrel solar pv calc for Yakima Wa %0.1973;% for pheonix Az % 0.1451 for Pullman Wa
InstalledSolar = (iSolar/CapacityFactor)/(365*24); %go from MWh to MW installed



RawData(:,1) = (iDemand/iWind)*RawData(:,1);%in MW increase wind to match demand
RawData(:,3) = (iDemand/iSolar)*RawData(:,3);
RawData(:,9) = (1-PercSolar)*RawData(:,1) + PercSolar*RawData(:,3);
%c = iWind/(365*24*4782) = 0.2741;

InstalledWind = (1-PercSolar)*(iDemand/iWind)*4782;%(in MW)  4782 from BPA website
InstalledSolar = PercSolar*(iDemand/iSolar)*InstalledSolar;
iWind = iDemand*(1-PercSolar);
iSolar = iDemand*PercSolar;
if isnan(InstalledWind)
    InstalledWind = 0;
end
if isnan(InstalledSolar)
    InstalledSolar = 0;
end


RawData(:,5) = RawData(:,2) - RawData(:,1);%net power as a function of time
RawData(:,6) = RawData(:,2) - RawData(:,3);% in MW
RawData(:,7) = RawData(:,2) - RawData(:,9);

iStorageBoth = (300/(60*60))*sum((RawData(:,7) > 0).*RawData(:,7));%amount of storage needed to supply demand in MWh

iLossesBoth = (iStorageBoth/StorageEff - iStorageBoth);%Amount of additional energy needed by storage efficiency losses
iFractionBoth = iLossesBoth/iDemand;

RawData(:,9) = RawData(:,9)*(1+iFractionBoth);
RawData(:,7) = RawData(:,2) - RawData(:,9);%net power as a function of time

InstalledSolar = (1+iFractionBoth)*InstalledSolar;
InstalledWind = (1+iFractionBoth)*InstalledWind;
%costs from NREL Energy Technology Cost and Performance Data


iWind = iWind*(1+iFractionBoth);
iSolar = iSolar*(1+iFractionBoth);
iStorageS = 0;
iStorageMax = 9e6;
%CostperMWhHydrogenTank = 14.69e3;%price per MWh, for hydrogen storage at 700 bar from DOE FY 2015 annual progress report 
for j = 1:100% converges towards proper energy balance
    iStorageBothT = 0*RawData(:,7);
    iDemandSatisfied1 = 0;
    for i = 1:length(RawData(:,5))% energy storage in MWh
        if i > 1
            if RawData(i,7) > 0
                iDemandSatisfied1 = iDemandSatisfied1 + 300/(60*60)*RawData(i,7);
                iStorageBothT(i) = iStorageBothT(i-1) - 300/(60*60)*RawData(i,7)/Ebat;
            else
                %iDemandSatisfied1 = iDemandSatisfied1;
                iStorageBothT(i) = iStorageBothT(i-1) - 300/(60*60)*RawData(i,7)*Ebat;
            end
        else
            if RawData(i,7) > 0
                iDemandSatisfied1 = 300/(60*60)*RawData(i,7);
                iStorageBothT(i) =  iStorageS -300/(60*60)*RawData(i,7)/Ebat;

            else
                iDemandSatisfied1 = 0;
                iStorageBothT(i) =  iStorageS -300/(60*60)*RawData(i,7)*Ebat;
            end
        end
        if j > 1
            %iStorageBothT(i) = iStorageBothT(i)*(1-((iStorageBothT(i)/iStorageMax))*(300*SelfDischarge));
        end
    end

    iStorageBothT = (iStorageBothT - min(iStorageBothT))/(iDemand);
    iStorageMax = max(iStorageBothT)*iDemand;
    iStorageS = iStorageBothT(1);
    EnergySurplus = (iStorageBothT(end) - iStorageBothT(1));
    RawData(:,9) = RawData(:,9)*(1-EnergySurplus*(2/3));
    RawData(:,7) = RawData(:,2) - RawData(:,9);
    InstalledSolar = InstalledSolar*(1-EnergySurplus*(2/3));
    InstalledWind = InstalledWind*(1-EnergySurplus*(2/3));
    iWind = iWind*(1-EnergySurplus*(2/3));
    iSolar = iSolar*(1-EnergySurplus*(2/3));
end

Analysis.BatterySize = iStorageMax;

n = 30;
DiscountRate = 0.056;% from eia cost determinations
DegredationRate = 0.01;%1% per year
OandM = 10;%assume roughly comparable to wind/solar power. $/MWh
TOandM = OandM*iDemand;

BatCost = Analysis.BatterySize*CostB;



for i = 1:n
    if i == 1
        LCOE1 = BatCost;
        LCOE2 = 0;
        LCOE2nodegrade = 0;
        LCOEH = 0;
    end
    LCOE1 = LCOE1 ;%+ TOandM/(1+DiscountRate)^i;
    LCOE2 = LCOE2 + (iDemand*(1-DegredationRate)^i)/((1+DiscountRate)^i);
    LCOE2nodegrade = LCOE2nodegrade + iDemand/((1+DiscountRate)^i);
end
LCOEbat = LCOE1/LCOE2;

Analysis.Batcost = LCOEbat;%this is a cost estimate because there is no lcoe data

CostperMWhSolarUnSub = 74.2;% from eia 2015 energy generation pdf for estimated lcoe for 2022
CostperMWhSolarSub = 58.2;%in $/MWh
altCostSolar = (CostperMWhSolarSub*iSolar/iDemand);
% CostperMWsolarUnSub = 0.26*24*365*CostperMWhSolarUnSub;%in $/MW installed
% CostperMWsolarSub = 0.26*24*365*CostperMWhSolarSub;
% SolarCostSub = InstalledSolar*CostperMWsolarSub/iDemand;%in $/MWh
% SolarCost = InstalledSolar*CostperMWsolarUnSub/iDemand;
Analysis.SolarCost = altCostSolar;


CostperMWhWindUnSub = 58.5;
CostperMWhWindSub = 50.9;
% CostperMWWindUnSub = 0.42*24*365*CostperMWhWindUnSub;
% CostperMWWindSub = 0.42*24*365*CostperMWhWindSub;
% WindCostSub = InstalledWind*CostperMWWindSub/iDemand;
% WindCost = InstalledWind*CostperMWWindUnSub/iDemand;

altCostWind = (CostperMWhWindSub*iWind/iDemand);


Analysis.WindCost = altCostWind;
end
