function [  ] = Figure6_4(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

A = 0.5;%percentage of demand met by baseline suppy. max value = 0.707
imax = 150;
jmax = 150;
maxScale = 3.6;

j = 1;

PercSolar = ((j-1)/(jmax-1));
% my method
for i = 1:imax
    Scale(i) = maxScale*(i/imax);
    if i == 1
        BatterySizing = 0;
        dBattery = 100;
    else
        BatterySizing = 0;
        dBattery = BatterySizingMatrix(j,i-1) + 100;
    end

    [LCOEtemp,Penetrationtemp,~,~] = RenewablesWithBattery(PercSolar,Scale(i),A,BatterySizing);
    check = 0;
    tol = 1e-6;
    while check == 0
        BatterySizing = BatterySizing + dBattery;
        [LCOEtemp2,Penetrationtemp2,~,~] = RenewablesWithBattery(PercSolar,Scale(i),A,BatterySizing);
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
    [LCOEtemp,PenetrationTemp,WindTemp,SolarTemp] = RenewablesWithBattery(PercSolar,Scale(i),A,BatterySizing);
    LCOEMatrix(j,i) = LCOEtemp;
    PenetrationMatrix(j,i) = PenetrationTemp;
    iWindMatrix(j,i) = WindTemp;
    iSolarMatrix(j,i) = SolarTemp;
end

%nadias dumb method
j = 2;
for i = 1:imax
    Scale(i) = maxScale*(i/imax);
    if i == 1
        BatterySizing = 0;
        dBattery = 100;
    else
        BatterySizing = BatterySizingMatrix(j,i-1);
        dBattery = BatterySizingMatrix(j,i-1)/10 + 100;
    end

    [LCOEtemp,Penetrationtemp,~,~] = RenewablesWithBattery(PercSolar,Scale(i),A,BatterySizing);
    check = 0;
    tol = 1e-6;
    while check == 0
        BatterySizing = BatterySizing + dBattery;
        [LCOEtemp2,Penetrationtemp2,~,~] = RenewablesWithBattery(PercSolar,Scale(i),A,BatterySizing);
        if (LCOEtemp2-LCOEtemp)/(Penetrationtemp2-Penetrationtemp) > 0
            BatterySizing = BatterySizing - dBattery; 
            dBattery = dBattery/5;
            if abs(LCOEtemp2-LCOEtemp)/LCOEtemp < tol
                BatterySizing = BatterySizing + 5*dBattery;
                check = 1;
                disp(i);
            end
        elseif (LCOEtemp2-LCOEtemp)/(Penetrationtemp2-Penetrationtemp) < 0
            LCOEtemp = LCOEtemp2;
            Penetrationtemp = Penetrationtemp2;
        end
    end
    BatterySizingMatrix(j,i) = BatterySizing;
    [LCOEtemp,PenetrationTemp,WindTemp,SolarTemp] = RenewablesWithBattery(PercSolar,Scale(i),A,BatterySizing);
    LCOEMatrix(j,i) = LCOEtemp;
    PenetrationMatrix(j,i) = PenetrationTemp;
    iWindMatrix(j,i) = WindTemp;
    iSolarMatrix(j,i) = SolarTemp;
end

j = 3;
%no energy storage
for i = 1:imax
    Scale(i) = maxScale*(i/imax);
    BatterySizingMatrix(j,i) = 0;
    [LCOEtemp,PenetrationTemp,~,~,~] = CurtailedRenewables(PercSolar,Scale(i),A);
    LCOEMatrix(j,i) = LCOEtemp;
    PenetrationMatrix(j,i) = PenetrationTemp;
end




figure(1)
plot(PenetrationMatrix(1,:),LCOEMatrix(1,:),'LineWidth',2);
hold on
plot(PenetrationMatrix(2,:),LCOEMatrix(2,:),'LineWidth',2);
plot(PenetrationMatrix(3,:),LCOEMatrix(3,:),'LineWidth',2);
ylabel('LCOE ($/MWh)')
xlabel('Annual Renewable Energy Penetration');
legend('Penetration Maximization','LCOE Minimization','No Energy Storage');
end
