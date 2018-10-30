function [ Flow, Hydrogen, Eff ] = HTM2( Flow,PInlet,PHyd )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
PR = PInlet/PHyd;

H2Selectivity = 5000;%H2/N2  assume that it's the same ratio for everything
A = 0.5;
NetFlow = Flow.H2 + Flow.CH4 + Flow.CO + Flow.CO2 + Flow.H2O;
X = Flow.H2/NetFlow;

%assume pressure doesn't change along length of HTM

PPInlet = X*PInlet;
PPHyd = PHyd*H2Selectivity/(H2Selectivity+1);

if PPHyd >= PPInlet
    Hydrogen.H2 = 1e-30*Flow.H2;
    Hydrogen.CO = 0;
    Hydrogen.CO2 = 0;
    Hydrogen.H2O = 0;
    Hydrogen.CH4 = 0;
    Hydrogen.T = Flow.T;
    Eff = 0;
else
    PPOutlet = PPHyd;
    Hydrogen.H2 = A*(Flow.H2*(PInlet/PPHyd)-NetFlow)/((PInlet/PPHyd)-1);
    
    NetFlow = Flow.H2O + Flow.CH4 + Flow.CO + Flow.CO2;
    
    Flow.H2 = Flow.H2 - Hydrogen.H2;
    Hydrogen.CO2 = Hydrogen.H2/H2Selectivity*Flow.CO2/NetFlow;
    Flow.CO2 = Flow.CO2 - Hydrogen.CO2;
    Hydrogen.CO = Hydrogen.H2/H2Selectivity*Flow.CO/NetFlow;
    Flow.CO = Flow.CO - Hydrogen.CO;
    Hydrogen.CH4 = Hydrogen.H2/H2Selectivity*Flow.CH4/NetFlow;
    Flow.CH4 = Flow.CH4 - Hydrogen.CH4;
    Hydrogen.H2O = Hydrogen.H2/H2Selectivity*Flow.H2O/NetFlow;
    Flow.H2O = Flow.H2O - Hydrogen.H2O;
    
    Hydrogen.T = Flow.T;
    Eff = Hydrogen.H2/(Flow.H2 + Hydrogen.H2);
end

end

