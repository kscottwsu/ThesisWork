function [RenewablesProfile,CurtailableDemand,iSolar,iWind,MaxRamp] = ProfileGeneration(PercSolar,A)
RawWindData = load('WindData.mat');
RawData = RawWindData.RawWindData.Spring;
RawData(105121:end,:) = [];%trim length to one year


RawSolarData = load('SolarData.mat');
RawSolarData = RawSolarData.SolarData;


n = 0:(length(RawData(:,1))-1);
t = n*5*60;%convert from 5 minutes to seconds
t = t/(60*60*24);%in days
RawData(:,4) = t';%time in seconds

SolarPowerGraphScript;

Demand = (5*60/(60*60))*RawData(:,2);%in MWh
iDemand = sum(Demand);

iWind = 300/(60*60)*sum(RawData(:,1));%roughly 4782 MW installed capacity
iSolar = 300/(60*60)*sum(RawSolarData);%solar data in MWh
    
WindData = 300/(60*60)*RawData(:,1)*(iDemand)/iWind;
SolarData = 300/(60*60)*RawSolarData*(iDemand)/iSolar;



RenewablesProfile = PercSolar*SolarData + (1-PercSolar)*WindData;
iSolar = PercSolar*iDemand;
iWind = (1-PercSolar)*iDemand;
   

BaselineGeneration = A*sum(Demand)/length(Demand);%in MWh
BaselineGeneration = BaselineGeneration*ones(length(Demand),1);

CurtailableDemand = Demand - BaselineGeneration;
MaxTRampRate = max(CurtailableDemand)*ThermalLimitCalculation();

CDRR = (CurtailableDemand(2:end) - CurtailableDemand(1:(end-1)))/max(5*CurtailableDemand);
MaxDRampRate = max(max(CDRR),-min(CDRR));
%MaxRamp = MaxDRampRate*max(CurtailableDemand);

MaxRamp = MaxTRampRate;
% RMaxamp = 0.02*max(CurtailableDemand);
% 
% 
% figure(1)
% plot(t,((60^2)/300)*Demand);
% hold on
% plot(t,((60^2)/300)*BaselineGeneration,'LineWidth',2)
% hold off
% xlabel('Time (Days)');
% ylabel('MW');
% ylim([0 ((60^2)/300)*1.1*max(Demand)]);
% xlim([0 max(t)]);
% legend('Demand','Baseline Generation');

end

