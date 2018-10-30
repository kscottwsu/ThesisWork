function [Out1,Out2] = enthalpy(varargin) % enthalpy (h) and sensible enthalpy(h_s), if species and flow are provided it returns H (total enthalpy(kW)) and H_s Sensible enthalpy (kW).
% returns both total and sensible enthalpy in kJ/kmol, or in rate form kJ/s
% Option 1: provide a vector of temperatures and it returns total and specific enthalpy at those temperatures for all species CH4, CO, CO2, H2, H2O, N2, O2, C, NO, OH, H
% Option 2: provide a vector of temperatures and a cell array of strings for the species of interest
% Option 3: provide a structure where __.T coresponds to temperature, ___.CH4 coresponds to the flow rate of methane ____.H2 to the flow rate of hydrogen...
% Option 3, returns the rate of energy flow in (kJ/s)
if isstruct(varargin{1})
    Inlet=  varargin{1};
    T = Inlet.T;
    spec = fieldnames(Inlet);
    spec = spec(~strcmp('T',spec)); 
elseif isnumeric(varargin{1})
    T = varargin{1};
    if length(varargin)==2
        spec = varargin{2};
    else
        spec = {'CH4';'CO';'CO2';'H2';'H2O';'N2';'O2';'H';'OH';'C';'NO';};
    end
end

T1 = (T/1000);
T2 = ((T/1000).^2)/2; 
T3 = ((T/1000).^3)/3; 
T4 = ((T/1000).^4)/4;
T5 = -1./(T/1000); 
T6 = ones(length(T),1);

for i = 1:1:length(spec)
    h298 = 0;
    switch spec{i}
        case 'CH4'
            A = [85.81217,11.26467,-2.114146,0.138190,-26.42221,-153.5327,224.4143];
            B = [-0.703029,108.4773,-42.52157,5.862788,0.678565,-76.84376,158.7163];
            C = (T>1300)*A+(T<=1300)*B;
            hf = -74873.1;
        case 'CO'
            A = [35.15070,1.300095,-.205921,0.013550,-3.282780,-127.8375,231.7120];
            B = [25.56759,6.096130,4.054656,-2.671301,0.131021,-118.0089,227.3665;];
            C = (T>1300)*A+(T<=1300)*B;
            hf = -110527.1;
        case 'CO2'
            A = [58.16639,2.720074,-0.492289,0.038844,-6.447293,-425.9186,263.6125;];
            B = [24.99735,55.18696,-33.69137,7.948387,-0.136638,-403.6075,228.2431;];
            C = (T>1200)*A+(T<=1200)*B;
            hf = -393522.4;
        case 'H2'
            A = [43.413560,-4.293079,1.272428,-.096876,-20.533862,-38.515158,162.081354;];
            B = [18.563083,12.257357,-2.859786,0.268238,1.977990,-1.147438,156.288133;];
            D = [33.066178,-11.363417,11.432816,-2.772874,-0.158558,-9.980797,172.707974;];
            C = (T>2500)*A+(T<=2500).*(T>1000)*B+(T<=1000)*D;
            hf = 0;
        case 'H2O'
            A = [41.96426,8.622053,-1.499780,0.098119,-11.15764,-272.1797,219.7809;];
            B = [30.09200,6.832514,6.793435,-2.534480,0.082139,-250.8810,223.3967;];
            C = (T>1700)*A+(T<=1700)*B;
            hf = -241826.4; %water vapor
        case 'N2'
            A = [35.51872,1.128728,-0.196103,0.014662,-4.553760,-18.97091,224.9810;];
            B = [19.50583,19.88705,-8.598535,1.369784,0.527601,-4.935202,212.3900;];
            D = [28.98641,1.853978,-9.647459,16.63537,0.000117,-8.671914,226.4168;];
            C = (T>2000)*A+(T<=2000).*(T>500)*B+(T<=500)*D;
            hf = 0;
        case 'O2'
            A = [20.91111,10.72071,-2.020498,0.146449,9.245722,5.337651,237.6185;];
            B = [30.03235,8.772972,-3.988133,0.788313,-0.741599,-11.32468,236.1663;];
            D = [31.32234,-20.23531,57.86644,-36.50624,-0.007374,-8.903471,246.7945;];
            C = (T>2000)*A+(T<=2000).*(T>700)*B+(T<=700)*D;
            hf = 0;
        case 'C'
            C = [21.1751,-0.812428,0.448537,-0.043256,-0.013103,710.347,183.8734];
            hf = 716669;
        case 'NO'
            A = [35.99169,0.95717,-0.148032,0.009974,-3.004088,73.10787,246.1619];
            B = [23.83491,12.58878,-1.139011,-1.497459,0.21419,83.35783,237.1219];
            C = (T>1200)*A+(T<=1200)*B;
            hf = 90291.14;
        case 'OH'
            A = [28.74701,4.7144,-0.814725,0.054748,-2.747829,26.41439,214.1166];
            B = [32.27768,-11.36291,13.60545,-3.846486,-0.001335,29.75113,225.5783];
            C = (T>1300)*A+(T<=1300)*B;
            hf = 38987.06;
        case 'H'
            C =[20.78603,4.85E-10,-1.58E-10,1.53E-11,3.20E-11,2.12E+02,1.40E+02];
            hf = 218194;
        case 'C2H6'
            C =[6.9,172.7,-64.06,7.285,9.173,0,0,0];
            h298 = 9173.70215 + 83700;
            hf = -83700;
        case 'C3H8'
            C =[-4.04,304.8,-157.2,31.74,11.05,0,0,0];
            h298 = 11005.6940 +104700;
            hf = -104700;
        case 'C6H6'
            C =[-50.24,568.2244,-442.503,134.5489,6.6206,0,0,0];
            h298 = 6620.6417 + 82930;
            hf = -82930;%benzene vapor
    end
    h.(spec{i}) = (T1.*C(:,1) + T2.*C(:,2) + T3.*C(:,3) + T4.*C(:,4) + T5.*C(:,5) + T6.*C(:,6))*1000-h298;%convert from kJ/mol to kJ/kmol
    h_s.(spec{i}) = h.(spec{i}) - hf;
end

if ~exist('Inlet','var')
    Out1 = h;
    Out2 = h_s;
else
    Out1 = 0;
    Out2 = 0;
    for i = 1:1:length(spec)
        Out1 = Out1 + h.(spec{i}).*Inlet.(spec{i});
        Out2 = Out2 + h_s.(spec{i}).*Inlet.(spec{i});
    end
end