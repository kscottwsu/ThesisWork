function [ ] = reSOFCPlottingCell( )

SOFCModePlotting(1,101,1);
SOFCModePlotting(1,2000,2);
SOFCModePlotting(0,2000,3);

figure(10)
ylabel('Voltage (V)');
xlabel('Current Density (A/cm^2)');
legend('Unpressurized H2','Pressurized H2','Pressurized CH4');

figure(11)
ylabel('Power Density (W/cm^2)');
xlabel('Current Density (A/cm^2)');
legend('Unpressurized H2','Pressurized H2','Pressurized CH4');

figure(12)
ylabel('Efficiency')
xlabel('Current Density (A/cm^2)');
legend('Unpressurized H2','Pressurized H2','Pressurized CH4','Oxidation-reduction Efficiency');

SOECModePlotting(1,101,4);
SOECModePlotting(1,2000,5);
SOECModePlotting(0,2000,6);

figure(13)
ylabel('Voltage (V)');
xlabel('Current Density (A/cm^2)');
legend('Unpressurized H2','Pressurized H2','Pressurized CH4');

figure(14)
ylabel('Power Density (W/cm^2)');
xlabel('Current Density (A/cm^2)');
legend('Unpressurized H2','Pressurized H2','Pressurized CH4');

figure(15)
ylabel('Efficiency')
xlabel('Current Density (A/cm^2)');
legend('Unpressurized H2','Pressurized H2','Pressurized CH4','Oxidation-reduction Efficiency');
end