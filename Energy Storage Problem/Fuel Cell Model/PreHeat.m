function [ Flow, WPreHeat ] = PreHeat( varargin )


if length(varargin) == 3
    Flow = varargin{1};
    Water = varargin{2};
    Tfc = varargin{3};

    [~,Hw1] = enthalpy(Water);
    Water.T = Tfc;
    [~,Hw2] = enthalpy(Water);
    [~,Hf1] = enthalpy(Flow);
    Flow.T = Tfc;
    [~,Hf2] = enthalpy(Flow);

    WPreHeat = ((Hw2 - Hw1) + (Hf2 - Hf1))/1000;%in KW

    Flow.H2O = Water.H2O;

else
    Flow = varargin{1};
    Tfc = varargin{2};
    
    [~,H1] = enthalpy(Flow);
    Flow.T = Tfc;
    [~,H2] = enthalpy(Flow);
    WPreHeat = (H2-H1)/1000;
    
end


end