inlet = 1;
utilization = 0.7;
recirculation = 0.75;
anodeStream = 0;
for i = 1:1e9
    anodeStream = anodeStream+inlet;
    anodeStream = anodeStream*(1-utilization);
    exhaust(i) = anodeStream*(1-recirculation);
    anodeStream = anodeStream*recirculation;
end
plot(exhaust);
disp(exhaust(end));