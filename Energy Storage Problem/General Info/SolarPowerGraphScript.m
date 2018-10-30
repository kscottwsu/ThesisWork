%% Solar Day Plotting Script

n1 = find(t == 127);
n2 = find(t == 136);

%data = RawSolarData(n1:n2);
data = RawData(n1:n2,1);
figure(1);
plot(t(n1:n2),data);
xlabel('Time (Day of Year)');
ylabel('Power Output (MW)');


