function HTMOptimizationPlot( ~ )
%shows HTM pressure vs system efficiency and hydrogen recovery mol flow
%rate

%calculated at 0.695 Amp/cm^2


FlowIn = 1e-6;

PHTMmin = 1;
PHTMmax = 150;

imax = 100;


for i = 1:imax
    tic
    PHTM(i) = PHTMmin + (PHTMmax - PHTMmin)*(i-1)/(imax-1);
    [~,ElecEff(i),HydProduction(i),j(i),Vfc(i),FCEff(i),PFuelCell(i),NetEff(i),~] = FuelCellSystemModel( FlowIn,3,PHTM(i),2000,0,0.7,0);
    toc
end


disp(PHTM(NetEff == max(NetEff)));%display pressure at which net efficiency is maximized

figure(1)
yyaxis left
plot(PHTM,ElecEff,'LineWidth',2);
hold on
plot(PHTM,NetEff,'LineWidth',2);
hold off
%ylim([0.6 0.8])
xlabel('HTM Pressure (KPa)');
ylabel(' System Efficiency');
yyaxis right
plot(PHTM,HydProduction,'LineWidth',2);
ylabel('H2 Production (Mols/S cm^2)');
hold off
xlim([PHTMmin PHTMmax])
legend('Electrical Efficiency','Net Efficiency','Hydrogen Production')

%net efficiency increases at low pressures due to higher amount of hydrogen
%recovery.  Electrical efficiency decreases at high pressures due to
%increase in heating demand and CO2 compressor demand
end

