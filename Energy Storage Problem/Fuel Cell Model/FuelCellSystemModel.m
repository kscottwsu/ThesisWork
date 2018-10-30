function [EOut,SysEff,HydProduction,j,Vfc,FCEff,PFuelCell,NetSysEff,OxiEff] = FuelCellSystemModel( FlowIn,DS2CRatio,PHTM,PFC,bHyd,FCUtilization,bDisp)
Tfc = 923.15;%650 C
F = 96485.332;%faraday constant C/mol

%FlowIn = mol/sec CH4
%DFuelUtilization = desired fuel utilization rate, used to control fuel
%recycling valve.  Needs to be above 70%
%DS2CRatio = desired steam to carbon ratio, used to ensure that carbon
%isn't deposited on the fuel cell

%Only has ohmic losses, so results unrealistic past j = ~1 A/cm^2

h = enthalpy(298);
HMethane = -(2*h.H2O + h.CO2 -h.CH4-2*h.O2); %2H2 + O2 -->  2H2O, ion transfer
HHyd = -(h.H2O-h.H2-h.O2/2);
if bHyd ~= 1
    Flow.CH4 = FlowIn;% in mols/sec
    Flow.CO = 0;
    Flow.CO2 = 0;
    Flow.H2 = 0;
    Flow.H2O = 0;
    Flow.T = 333.15;%in K
    EIn = HMethane*Flow.CH4/1000;
else
    Flow.CH4 = 0;
    Flow.CO = 0;
    Flow.CO2 = 0;
    Flow. H2 = 4*FlowIn;
    Flow.H2O = 0;
    Flow.T = 333.15;
    EIn = HHyd*Flow.H2/1000;
end
    
    FlowPreHeat = Flow;
    [Flow, ~] = PreHeat(Flow, Tfc);
    Inlet = Flow;
    RecircPerc = 0.75;
    imax = 100;
    for i = 1:imax
        if i == 1
            Flow = Inlet;
            Water.H2O = Flow.CH4*DS2CRatio;
            Water.T = Flow.T;
            Flow = FlowAdd(Flow,Water);
        elseif i~=1
            Flow = FlowAdd(Flow,Inlet);
        end

        if i~=imax
            rxn = FCUtilization*(Flow.CO + Flow.H2 + 4*Flow.CH4);
            [Flow,ErrorFlag] = FCellLite(Flow,rxn);
        else
            [Flow,Oxy, PFuelCell, QFuelCell, Vfc,j,FCEff,OxiEff,ErrorFlag] = FuelCell(Flow,PFC,FCUtilization,bDisp);
        end
        if ErrorFlag == 1
            PFuelCell = 0;
        end
        Oxy.T = 333.15;
        temp = Flow;
        Outlet = FlowMult(temp,(1-RecircPerc));
        Flow = FlowMult(temp,RecircPerc);
    end
    if bHyd ~= 1
        FuelUtilization = (j/(2*F))/(4*Inlet.CH4);
    else
        FuelUtilization = (j/(2*F))/(Inlet.H2);
    end
    PFuelCell = PFuelCell/1000;
    CS2CRatio = Flow.H2O/(Flow.CO/2 + Flow.CH4 + Inlet.CH4);
    errorS2C = (DS2CRatio - CS2CRatio);%Amount of water needed relative to the amount of carbon flowing in to ensure SC2 Ratio is good
    errorS2C = max(0,errorS2C);
    Water.H2O = errorS2C*FlowIn;%amount of additional water needed to maintain S2C Ratio
    Water.T = 298;% assume water starts at STP
    [Water, QEvaporator] = Evaporator(Water);%boil water to 100 C, doesn't work this way...
    [Water, PCompressor] = Compressor(Water, 101, PFC);%compress up to 20 bar
    if isnan(PCompressor)
        PCompressor = 0;
        Water.T = 0;
    end
    WaterComp = Water;
    
    [Outlet,Hydrogen,Eff] = HTM2(Outlet,PFC,PHTM);%Use high temperature hydrogen transport membrane to extract hydrogen fuel  
    Hydrogen1 = Hydrogen;%heat exchanger with fuel
    [Hydrogen,~] = HeatExchanger(Hydrogen,373);%Heat exchanger with water to remove remainder of heat after going through heat exchanger
    
    PR = (700*100/PHTM)^(1/3);
    [Hydrogen,PHydCompressor1] = Compressor(Hydrogen,PHTM,PHTM*PR);
    
    [Hydrogen,~] = HeatExchanger(Hydrogen,373);%Heat exchanger with water to remove remainder of heat after going through heat exchanger
    [Hydrogen,PHydCompressor2] = Compressor(Hydrogen,PHTM,PHTM*PR);
    
    [Hydrogen,~] = HeatExchanger(Hydrogen,373);%Heat exchanger with water to remove remainder of heat after going through heat exchanger
    [Hydrogen,PHydCompressor3] = Compressor(Hydrogen,PHTM,PHTM*PR);
    PHydCompressor = PHydCompressor1+PHydCompressor2+PHydCompressor3;
    Hydrogen2 = Hydrogen;%heat exchanger with water, then heat exchanger with oxygen inlet
    Hydrogen.T = 298;%heat exchanger with air
    [Hydrogen,QWGSHyd] = WGSReactor(Hydrogen);
    
    TMax = 1200;%max turbine intake temperature in K
    [Outlet,OxyAdd] = Oxidizer(Outlet,TMax);
    if Outlet.T > TMax
        Outlet.T = TMax;
    end
    
    [Outlet,PTurbine] = Turbine(Outlet,PFC,101);
    [Outlet,QHXflow] = HeatExchanger(Outlet,373);%Used to boil water in Evaporator
    [Outlet,QCondenser] = Condenser(Outlet,298);%Used to boil water in Evaporator
    if bHyd ~= 1
        [Outlet, QWGSCO2] = WGSReactor(Outlet);
        PR = (700*100/101)^(1/3);
        [Outlet, PCO2Compressor1] = Compressor(Outlet, 101, 101*PR);
        [Outlet,~] = HeatExchanger(Outlet,373);%Heat exchanger with water to remove remainder of heat after going through heat exchanger

        [Outlet, PCO2Compressor2] = Compressor(Outlet, 101, 101*PR);
        [Outlet,~] = HeatExchanger(Outlet,373);%Heat exchanger with water to remove remainder of heat after going through heat exchanger

        [Outlet, PCO2Compressor3] = Compressor(Outlet, 101, 101*PR);
        [Outlet,~] = HeatExchanger(Outlet,373);%Heat exchanger with water to remove remainder of heat after going through heat exchanger
        
        PCO2Compressor = PCO2Compressor1+PCO2Compressor2+PCO2Compressor3;
        
        CO2 = Outlet;%CO2 stream is run through heat exchanger with Fuel Flow
        Outlet.T = 298;%heat exchanger with air 
        [FlowPreHeata, ~] = HeatExchanger2(FlowPreHeat,CO2);

        QEvaporation = max(0,QEvaporator - QCondenser - QHXflow);
        [~,QFuelHeat] = PreHeat(FlowPreHeata, Tfc);
        if QFuelHeat < 0
            QFuelHeat = 0;
        end
        [~,QWaterHeat] = PreHeat(WaterComp,Tfc);
        if isnan(QWaterHeat)
            QWaterHeat = 0;
        end

        Oxy.O2 = Oxy.O2 + OxyAdd.O2;

        [Oxya, ~] = HeatExchanger2(Oxy,Hydrogen1);%Hydrogen Stream is run though heat exchanger with Oxy Stream
        [Oxyb, ~] = HeatExchanger2(Oxya,Hydrogen2);%Compressed Hydrogen Stream is run through heat exchanger with Oxy Stream
        [~,QOxyHeat] = PreHeat(Oxyb,Tfc);
        PHeating = max(0,QFuelHeat + QWaterHeat + QOxyHeat + QEvaporation - QFuelCell)/1000;

        EOut = PTurbine + 0.98*PFuelCell - PCompressor - PCO2Compressor - PHydCompressor - PHeating;
    else
        [~,QFuelHeat] = PreHeat(FlowPreHeat,Tfc);
        Oxy.O2 = Oxy.O2 + OxyAdd.O2;

        [Oxya, ~] = HeatExchanger2(Oxy,Hydrogen1);%Hydrogen Stream is run though heat exchanger with Oxy Stream
        [Oxyb, ~] = HeatExchanger2(Oxya,Hydrogen2);%Compressed Hydrogen Stream is run through heat exchanger with Oxy Stream
        [~,QOxyHeat] = PreHeat(Oxyb,Tfc);
        PHeating = max(0,QFuelHeat  + QOxyHeat - QFuelCell)/1000;
        EOut = PTurbine + 0.98*PFuelCell - PCompressor - PHydCompressor - PHeating;
    end
    
    h = enthalpy(298);
    hrxn = -(h.H2O-h.H2-h.O2/2); %2H2 + O2 -->  2H2O, ion transfer
    hrxn2 = -(2*h.H2O + h.CO2 - h.CH4 - 2*h.O2);
    HydEnergy = Hydrogen.H2*hrxn/1000;
    
    
SysEff = EOut/EIn;
NetSysEff = (EOut + HydEnergy)/EIn;
HydProduction = Hydrogen.H2;%in Mols
end

function [Flow2] = FlowMult(Flow,R)
list = fieldnames(Flow);

for i = 1:length(list)
    if list{i} == 'T'
        Flow2.T = Flow.T;
    else
        Flow2.(list{i}) = Flow.(list{i})*R;
    end
end

end

function [Flow2] = FlowAdd(Flow,Flow2)
%adds two flows together, but doesn't change temperatures
list = fieldnames(Flow);

for i = 1:length(list)
    if list{i} ~= 'T'
       if isfield(Flow2,list{i})
           Flow2.(list{i}) = Flow2.(list{i}) + Flow.(list{i});
       else
           Flow2.(list{i}) = Flow.(list{i});
       end
    end
        
        
end
end

function [Flow,ErrorFlag] = FCellLite(Flow,Rrxn1)
ErrorFlag = 0;

RrxnVec = 0;
check = 0;
Rrxn1Guess = Rrxn1/15;
[Flow,~] = ChemicalEquilibrium(Flow,20);

while check == 0
    if (Rrxn1Guess > Flow.H2*0.75)
        Rrxn1Guess = Rrxn1Guess/2;
    elseif (Rrxn1Guess>(Rrxn1-sum(RrxnVec)))
        Rrxn1Guess = Rrxn1Guess/2;
    else
        Flow.H2O = Flow.H2O + Rrxn1Guess;
        Flow.H2 = Flow.H2 - Rrxn1Guess;
        [Flow,~] = ChemicalEquilibrium(Flow,20);
        RrxnVec(end+1) = Rrxn1Guess;
    end
    
    if sum(RrxnVec)>=((1-1e-4)*Rrxn1)
        check = 1;
    end
end

end

