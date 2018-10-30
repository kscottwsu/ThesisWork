function [ ] = DemandProfilePlotting( )
A = 0.5;%Assumed value of baseload generation relative to total yearly energy demand

RawWindData = load('WindData.mat');
RawData = RawWindData.RawWindData.Spring;

n = 0:(length(RawData(:,1))-1);
t = n*5*60;%convert from 5 minutes to seconds
t = t/(60*60*24);%in days
RawData(:,4) = t';%time in seconds

Demand = RawData(:,2);%in MW

BaselineGeneration = A*sum(Demand)/length(Demand);%in MW
BaselineGeneration = BaselineGeneration*ones(length(Demand),1);

figure(1)
plot(t,Demand);
hold on
plot(t,BaselineGeneration);
hold off
title('Demand Profile');
xlabel('Time (days)');
ylabel('Power (MW)');
legend('Demand','Baseload Generation');
xlim([0,370]);
ylim([0,10000]);
end

