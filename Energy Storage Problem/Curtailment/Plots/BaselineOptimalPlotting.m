function [ ] = BaselineOptimalPlotting( )
imax = 200;
jmax = 100;
maxScale = 2;

A = 0.4;
for j = 1:jmax
    PercSolar(j) = 0.4*((j-1)/(jmax-1));
    for i = 1:imax
        Scale(i) = maxScale*(i/imax);
        [LCOEtemp,PenetrationTemp,WindTemp,SolarTemp] = CurtailedRenewables(PercSolar(j),Scale(i),A);   
        LCOEMatrix(j,i) = LCOEtemp;
        PenetrationMatrix(j,i) = PenetrationTemp;
        iWindMatrix(j,i) = WindTemp;
        iSolarMatrix(j,i) = SolarTemp;
    end
end

check = 0;
count = 0;
while check == 0
    count = count+1;
    if count == 1
       LCOEOptimalA(count) = min(LCOEMatrix(:,1));
       n = find(LCOEMatrix(:,1) == LCOEOptimalA(count));
       PenetrationOptimalA(count) = PenetrationMatrix(n,1);
       iWindOptimalA(count) = iWindMatrix(n,1);
       iSolarOptimalA(count) = iSolarMatrix(n,1);
       PercSolarOptimalA(count) = PercSolar(n);
    else
        dLCOEdPen = (LCOEMatrix - LCOEOptimalA(count-1))./(PenetrationMatrix - PenetrationOptimalA(count-1));
        if max(dLCOEdPen(:)) > 0
            dLCOEdPen(PenetrationMatrix < PenetrationOptimalA(count-1)) = -max(dLCOEdPen(:));
        else
            dLCOEdPen(PenetrationMatrix < PenetrationOptimalA(count-1)) = -1;
        end
        if max(max(dLCOEdPen)) > 0
            dLCOEdPen(dLCOEdPen<0) = 10*max(dLCOEdPen(:));
            [M,I] = min(dLCOEdPen(:));
            [i_row,i_col] = ind2sub(size(dLCOEdPen),I);
            LCOEOptimalA(count) = LCOEMatrix(i_row,i_col);
            PenetrationOptimalA(count) = PenetrationMatrix(i_row,i_col);
            iWindOptimalA(count) = iWindMatrix(i_row,i_col);
            iSolarOptimalA(count) = iSolarMatrix(i_row,i_col);
            PercSolarOptimalA(count) = PercSolar(i_row);
            if i_col == imax
                check = 1;
            end
        else
            check = 1;
        end
    end
end


A = 0.5;
for j = 1:jmax
    PercSolar(j) = 0.3*((j-1)/(jmax-1));
    for i = 1:imax
        Scale(i) = maxScale*(i/imax);
        [LCOEtemp,PenetrationTemp,WindTemp,SolarTemp] = CurtailedRenewables(PercSolar(j),Scale(i),A);   
        LCOEMatrix(j,i) = LCOEtemp;
        PenetrationMatrix(j,i) = PenetrationTemp;
        iWindMatrix(j,i) = WindTemp;
        iSolarMatrix(j,i) = SolarTemp;
    end
end

check = 0;
count = 0;
while check == 0
    count = count+1;
    if count == 1
       LCOEOptimalB(count) = min(LCOEMatrix(:,1));
       n = find(LCOEMatrix(:,1) == LCOEOptimalB(count));
       PenetrationOptimalB(count) = PenetrationMatrix(n,1);
       iWindOptimalB(count) = iWindMatrix(n,1);
       iSolarOptimalB(count) = iSolarMatrix(n,1);
       PercSolarOptimalB(count) = PercSolar(n);
    else
        dLCOEdPen = (LCOEMatrix - LCOEOptimalB(count-1))./(PenetrationMatrix - PenetrationOptimalB(count-1));
        if max(dLCOEdPen(:)) > 0
            dLCOEdPen(PenetrationMatrix < PenetrationOptimalB(count-1)) = -max(dLCOEdPen(:));
        else
            dLCOEdPen(PenetrationMatrix < PenetrationOptimalB(count-1)) = -1;
        end
        if max(max(dLCOEdPen)) > 0
            dLCOEdPen(dLCOEdPen<0) = 10*max(dLCOEdPen(:));
            dLCOEdPen(PenetrationMatrix < PenetrationOptimalB(count-1)) = 10*max(dLCOEdPen(:));
            [M,I] = min(dLCOEdPen(:));
            [i_row,i_col] = ind2sub(size(dLCOEdPen),I);
            LCOEOptimalB(count) = LCOEMatrix(i_row,i_col);
            PenetrationOptimalB(count) = PenetrationMatrix(i_row,i_col);
            iWindOptimalB(count) = iWindMatrix(i_row,i_col);
            iSolarOptimalB(count) = iSolarMatrix(i_row,i_col);
            PercSolarOptimalB(count) = PercSolar(i_row);
            if i_col == imax
                check = 1;
            end
        else
            check = 1;
        end
    end
end


A = 0.6;
for j = 1:jmax
    PercSolar(j) = 0.3*((j-1)/(jmax-1));
    for i = 1:imax
        Scale(i) = maxScale*(i/imax);
        [LCOEtemp,PenetrationTemp,WindTemp,SolarTemp] = CurtailedRenewables(PercSolar(j),Scale(i),A);   
        LCOEMatrix(j,i) = LCOEtemp;
        PenetrationMatrix(j,i) = PenetrationTemp;
        iWindMatrix(j,i) = WindTemp;
        iSolarMatrix(j,i) = SolarTemp;
    end
end

check = 0;
count = 0;
while check == 0
    count = count+1;
    if count == 1
       LCOEOptimalC(count) = min(LCOEMatrix(:,1));
       n = find(LCOEMatrix(:,1) == LCOEOptimalC(count));
       PenetrationOptimalC(count) = PenetrationMatrix(n,1);
       iWindOptimalC(count) = iWindMatrix(n,1);
       iSolarOptimalC(count) = iSolarMatrix(n,1);
       PercSolarOptimalC(count) = PercSolar(n);
    else
        dLCOEdPen = (LCOEMatrix - LCOEOptimalC(count-1))./(PenetrationMatrix - PenetrationOptimalC(count-1));
        if max(dLCOEdPen(:)) > 0
            dLCOEdPen(PenetrationMatrix < PenetrationOptimalC(count-1)) = -max(dLCOEdPen(:));
        else
            dLCOEdPen(PenetrationMatrix < PenetrationOptimalC(count-1)) = -1;
        end

       if max(max(dLCOEdPen)) > 0
            dLCOEdPen(dLCOEdPen<0) = 10*max(dLCOEdPen(:));
            dLCOEdPen(PenetrationMatrix < PenetrationOptimalC(count-1)) = 10*max(dLCOEdPen(:));
            [M,I] = min(dLCOEdPen(:));
            [i_row,i_col] = ind2sub(size(dLCOEdPen),I);
            LCOEOptimalC(count) = LCOEMatrix(i_row,i_col);
            PenetrationOptimalC(count) = PenetrationMatrix(i_row,i_col);
            iWindOptimalC(count) = iWindMatrix(i_row,i_col);
            iSolarOptimalC(count) = iSolarMatrix(i_row,i_col);
            PercSolarOptimalC(count) = PercSolar(i_row);
            if i_col == imax
                check = 1;
            end 
       else
            check = 1;
       end
    end
end


figure(1)
plot(PenetrationOptimalA,LCOEOptimalA,'-x');
hold on 
plot(PenetrationOptimalB,LCOEOptimalB,'-x');
plot(PenetrationOptimalC,LCOEOptimalC,'-x');
hold off
title('Optimal LCOE of Renewables Generation With Varied Baseline Generation');
xlabel('Annual Renewable Energy Penetration');
ylabel('Net LCOE ($/MWh)');
legend('40% Baseline Generation','50% Baseline Generation','60% Baseline Generation');

figure(2)
plot(PenetrationOptimalA,PercSolarOptimalA,'-x');
hold on 
plot(PenetrationOptimalB,PercSolarOptimalB,'-x');
plot(PenetrationOptimalC,PercSolarOptimalC,'-x');
hold off
title('Optimal Percentage of Solar Energy With Varied Baseline Generation');
xlabel('Annual Renewable Energy Penetration');
ylabel('Optimal Percent of Renewable Energy from Solar');
legend('40% Baseline Generation','50% Baseline Generation','60% Baseline Generation');
end