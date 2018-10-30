function [ D ] = DiffusionCoefficient( Flow,P )
Ru = 8.3145;
Porosity = 0.26;
Tortuosity = 3;
Rpore = 1e-6;%in meters

%% atomic diffusion volumes
AtomicVolume.H = 1.98;
AtomicVolume.O = 5.48;
AtomicVolume.C = 16.5;

%Molecular Masses
AtomicMass.H = 1.00794;
AtomicMass.C = 12.011;
AtomicMass.O = 15.999;

% order: H2, H2O, CO, CO2, CH4
%find molar fraction of bulk material
NetFlow = Flow.H2 + Flow.H2O + Flow.CO + Flow.CO2 + Flow.CH4;

Xi(1) = Flow.H2/NetFlow;
Xi(2) = Flow.H2O/NetFlow;
Xi(3) = Flow.CO/NetFlow;
Xi(4) = Flow.CO2/NetFlow;
Xi(5) = Flow.CH4/NetFlow;

%Find diffusion Coefficients
for i = 1:5
    for j = 1:5
        if i == 1%H2
            Vi = 2*AtomicVolume.H;
            Mi = 2*AtomicMass.H;
        elseif i == 2%H2O
            Vi = 2*AtomicVolume.H + AtomicVolume.O;
            Mi = 2*AtomicMass.H + AtomicMass.O;
        elseif i == 3%CO
            Vi = AtomicVolume.C + AtomicVolume.O;
            Mi = AtomicMass.C + AtomicMass.O;
        elseif i == 4%CO2
            Vi = AtomicVolume.C + 2*AtomicVolume.O;
            Mi = AtomicMass.C + 2*AtomicMass.O;
        elseif i == 5%CH4
            Vi = AtomicVolume.C + 4*AtomicVolume.H;
            Mi = AtomicMass.C + 4*AtomicMass.H;
        end
        if j == 1%H2
            Vj = 2*AtomicVolume.H;
            Mj  = 2*AtomicMass.H;
        elseif j == 2%H2O
            Vj = 2*AtomicVolume.H + AtomicVolume.O;
            Mj = 2*AtomicMass.H + AtomicMass.O;
        elseif j == 3%CO
            Vj = AtomicVolume.C + AtomicVolume.O;
            Mj = AtomicMass.C + AtomicMass.O;
        elseif j == 4%CO2
            Vj = AtomicVolume.C + 2*AtomicVolume.O;
            Mj = AtomicMass.C + 2*AtomicMass.O;
        elseif j == 5%CH4
            Vj = AtomicVolume.C + 4*AtomicVolume.H;
            Mj = AtomicMass.C + 4*AtomicMass.H;
        end
        
        Mij(i,j) = 2*(1/Mi + 1/Mj)^-1;
        Mij(i,j) = Mij(i,j);
        
        Dij(i,j) = 1.43e-3*Flow.T^1.75/(P/100*Mij(i,j)^0.5*(Vi^(1/3)+Vj^(1/3))^2);%in cm^2/s
        Dij(i,j) = 10^-4*Dij(i,j);
    end
    %apply mixing rule

    Di2 = 0;
    for j = 1:5
       if i ~= j
           Di2 = Di2 + Xi(i)/Dij(i,j);
       end
    end
    
    Di(i) = (1-Xi(i))/(Di2);%m^2/sec
    
    %knudsen diffusion, important when pores are small
    DKi(i) = (2/3)*Rpore*((8*Ru*Flow.T)/(pi*Mi/1000))^0.5;%in m^2/sec
    
     
    %calculate effective diffusion
    DEffi(i) = ((Tortuosity/Porosity)*(DKi(i)^-1 + Di(i)^-1))^-1;
    
end

%% calculate diffusion rate of O2

MO2 = 2*AtomicMass.O;
MN2 = 2*14;
MO2N2 = 2*(1/MO2 + 1/MN2)^-1;
VO2 = 2*5.48;
VN2 = 2*5.69;

DO2N2 = 1.43e-3*Flow.T^1.75/(P/100*MO2N2^0.5*(VO2^(1/3)+VN2^(1/3))^2);%in cm^2/s
DO2N2 = DO2N2*10^-4;
DKO2 = (2/3)*Rpore*((8*Ru*Flow.T)/(pi*MO2/1000))^0.5;

DO2Eff = ((Tortuosity/Porosity)*(DKO2^-1 + DO2N2^-1))^-1;


D.H2Eff = DEffi(1);
D.H2 = Di(1);
D.H2OEff = DEffi(2);
D.H2O = Di(2);
D.O2 = DKO2;
D.O2Eff = DO2Eff;
end

