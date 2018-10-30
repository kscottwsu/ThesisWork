function [ ] = GetSolarData( ~ )
load('SolarRaw.mat');

list = fieldnames(Solar);
n = numel(list);

t = length(Solar.(list{1}).PowerMW);

SolarData = zeros(t,1);

for i = 1:t%time step
    for j = 1:n% plant number
        SolarData(i) = SolarData(i) + Solar.(list{j}).PowerMW(i);
    end
end

save('SolarData.mat','SolarData');
end

