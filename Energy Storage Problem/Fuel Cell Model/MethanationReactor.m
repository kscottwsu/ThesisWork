function [FlowOut,QNet,CO2in] = MethanationReactor(FlowOut)
%Methanation reactor, assuming low temperature and high pressure

h = enthalpy(298);
hrxn1 = h.H2O-h.H2-h.O2/2; %H2 + O2/2 -->  H2O, ion transfer
hrxn2 = h.CO2+h.H2-h.CO-h.H2O; %CO + H20 --> CO2 + H2, Water Gas Shift
hrxn3 = 3*h.H2+h.CO-h.CH4-h.H2O; %CH4+H2O --> CO + 3H2, Methane reforming


%WGS, all CO2 should form CO
Rrxn2 = FlowOut.CO2;
QWGS = hrxn2*Rrxn2;

FlowOut.CO2 = 0;
FlowOut.CO = FlowOut.CO + Rrxn2;
FlowOut.H2 = FlowOut.H2 - Rrxn2;
FlowOut.H2O = FlowOut.H2O + Rrxn2;

%Methanation, all CO2/CO should be reacted
Rrxn3 = FlowOut.CO;

FlowOut.H2 = FlowOut.H2 - 3*Rrxn3;
FlowOut.CO = FlowOut.CO - Rrxn3;
FlowOut.CH4 = FlowOut.CH4 + Rrxn3;
FlowOut.H2O = FlowOut.H2O + Rrxn3;


QMethanation = hrxn3*Rrxn3;
%add in more CO2 to turn as much H2 into CH4 as possible
Rrxn23 = FlowOut.H2/4;
CO2in.CO2 = Rrxn23;
CO2in.T = 333.15;

FlowOut.H2 = FlowOut.H2 - 4*Rrxn23;
FlowOut.H2O = FlowOut.H2 + 2*Rrxn23;
FlowOut.CH4 = FlowOut.CH4 + Rrxn23;

QWGS = QWGS + Rrxn23*hrxn2;
QMethanation = QMethanation + Rrxn23*hrxn3;

QNet = (QWGS + QMethanation)/1000;


end

