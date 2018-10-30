 function [ EUse,Outlet,SysEff,j,Vfc, ECEff,PElectrolysis,OxiEff ] = ElectrolysisCellSystemModel(FlowRate,SteamUtilization,PrEC,bHyd,bDisp)
Tfc = 923.15;%650 C
F = 96485.332;%faraday constant C/mol

if bHyd == 0
    %desired composition to maximize CH4 at 75% utilization
    %can't exceed 75% utilization due to Carbon Depostion Region
    FlowEquib.CO2 = 0.019986;
    FlowEquib.CO = 0.006524;
    FlowEquib.H2 = 0.28891;
    FlowEquib.H2O = 0.3323;
    FlowEquib.CH4 = 0.35228;
    FlowEquib.T = Tfc;
    
    H2Production = 4*FlowEquib.CH4 + FlowEquib.CO + FlowEquib.H2;
    H2OIn = FlowEquib.H2O+ 2*FlowEquib.CH4 + FlowEquib.H2;
    Util = H2Production/H2OIn;
    CO2temp = (FlowEquib.CO2+FlowEquib.CO+FlowEquib.CH4);
    CO2Fraction = CO2temp/H2OIn;
else
    FlowEquib.CO2 = 0;
    FlowEquib.CO = 0;
    FlowEquib.H2 = SteamUtilization;
    FlowEquib.H2O = (1-SteamUtilization);
    FlowEquib.CH4 = 0;
    FlowEquib.T = Tfc;
    CO2Fraction = 0;
end

WaterIn = FlowRate;
CO2In = CO2Fraction*FlowRate;

Inlet.H2O = WaterIn;
Inlet.CO2 = CO2In;
Inlet.CO = 0;
Inlet.H2 = 0;
Inlet.CH4 = 0;
Inlet.T = Tfc;


%check for carbon deposition region
H = 2*FlowEquib.H2 + 2*FlowEquib.H2O + 4*FlowEquib.CH4;
O = 1*FlowEquib.H2O + 1*FlowEquib.CO + 2*FlowEquib.CO2;
C = 1*FlowEquib.CH4 + 1*FlowEquib.CO + 1*FlowEquib.CO2;
temp = H + O + C;

H = H/temp;
O = O/temp;
C = C/temp;

RecircPerc = 0.75;
imax = 100;
for i = 1:imax
    if i == 1
        Flow = Inlet;
    else
        Flow = FlowAdd(Flow,Inlet);
    end
    if i ~= imax
        Rrxn1 = SteamUtilization*Flow.H2O;
       [Flow,ErrorFlag] = ECellLite(Flow,Rrxn1);
    else
        [Flow,OxyOut,PElectrolysis, QElectrolysis, Vfc, j,ErrorFlag,ECEff,OxiEff] = ElectrolysisCell(Flow,PrEC,SteamUtilization,bDisp);
    end
    Outlet = FlowMult(Flow,(1-RecircPerc));
    Flow = FlowMult(Flow,RecircPerc);
    if ErrorFlag == 1
        PElectrolysis = 10000000*PElectrolysis;
    end
end


Water.H2O = Inlet.H2O;
Water.T = 298;
[Water, QEvaporator] = Evaporator(Water);
[Water, PWCompressor] = Compressor(Water, 101, PrEC);
[Water, QPreHeatW] = PreHeat(Water, Tfc);

[OxyOut,QOxyCooling] = PreHeat(OxyOut,298);
QOxyCooling = -QOxyCooling;
[OxyOut, POxyCompressor] = Compressor(OxyOut,PrEC,160*100);

CO2.CO2 = Inlet.CO2;
CO2.T = 333.15;
[CO2, QPreHeatC] = PreHeat(CO2,Tfc);

[Outlet, PTurbine] = Turbine(Outlet,PrEC,101);
[Outlet, QCooling] = PreHeat(Outlet,398);
QCooling = -QCooling;
if bHyd == 0
    [Outlet,QMethanation,CO2in] = MethanationReactor(Outlet);
    CO2in.CO2 = CO2in.CO2 + CO2.CO2;%total amount of CO2 used
end
if bHyd == 1
    CO2in.CO2 = 0;
    QMethanation = 0;
end
[Outlet, QCondenser] = Condenser(Outlet,298);
if bHyd ~= 1
    PR = (160*100/101)^(1/3);
    [Outlet,~] = HeatExchanger(Outlet,298);%Heat exchanger with air to remove remainder of heat after going through heat exchanger
    [Outlet,PCH4Comp1] = Compressor(Outlet,101,PR*100);
    
    [Outlet,~] = HeatExchanger(Outlet,298);%Heat exchanger with air to remove remainder of heat after going through heat exchanger
    [Outlet,PCH4Comp2] = Compressor(Outlet,101,PR*100);
    
    [Outlet,~] = HeatExchanger(Outlet,298);%Heat exchanger with air to remove remainder of heat after going through heat exchanger
    [Outlet,PCH4Comp3] = Compressor(Outlet,101,PR*100);
    
    PCH4Compressor = PCH4Comp1+PCH4Comp2+PCH4Comp3;
else
    PR = (700*100/101)^(1/3);
    [Outlet,~] = HeatExchanger(Outlet,298);%Heat exchanger with air to remove remainder of heat after going through heat exchanger
    [Outlet,PHydCompressor1] = Compressor(Outlet,101,100*PR);

    [Outlet,~] = HeatExchanger(Outlet,298);%Heat exchanger with air to remove remainder of heat after going through heat exchanger
    [Outlet,PHydCompressor2] = Compressor(Outlet,101,100*PR);

    [Outlet,~] = HeatExchanger(Outlet,298);%Heat exchanger with air to remove remainder of heat after going through heat exchanger
    [Outlet,PHydCompressor3] = Compressor(Outlet,101,100*PR);
    PCH4Compressor = PHydCompressor1+PHydCompressor2+PHydCompressor3;
end
[Outlet,QCH4HX] = PreHeat(Outlet,298);

QBoiling = max(0,QEvaporator - QCondenser - QMethanation - QOxyCooling + QCH4HX - QCooling);

PHeating = max(0,QBoiling  + QPreHeatW + QPreHeatC - QElectrolysis);

h = enthalpy(298);
HMethane = -(2*h.H2O + h.CO2 - h.CH4 - 2*h.O2);
HH2 = -(h.H2O - h.H2 - 0.5*h.O2);
%HMethane = 4*HH2;
EIn = (HMethane*Outlet.CH4 + HH2*Outlet.H2)/1000;
EOut = PElectrolysis/0.98 + PWCompressor + POxyCompressor + PCH4Compressor + PHeating - PTurbine;
EUse = EOut;
SysEff = EIn/EOut;
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

function [Flow,ErrorFlag] = ECellLite(Flow,Rrxn1)
ErrorFlag = 0;RrxnVec = 0;
check = 0;
Rrxn1Guess = Rrxn1/15;
[Flow,~] = ChemicalEquilibrium(Flow,20);

while check == 0
    if (Rrxn1Guess > Flow.H2O*0.75)
        Rrxn1Guess = Rrxn1Guess/2;
    elseif (Rrxn1Guess>(Rrxn1-sum(RrxnVec)))
        Rrxn1Guess = Rrxn1Guess/2;
    else
        Flow.H2O = Flow.H2O - Rrxn1Guess;
        Flow.H2 = Flow.H2 + Rrxn1Guess;
        [Flow,~] = ChemicalEquilibrium(Flow,20);
        RrxnVec(end+1) = Rrxn1Guess;
    end
    
    if sum(RrxnVec)>=((1-1e-4)*Rrxn1)
        check = 1;
    end
end

end

