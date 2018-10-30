function [ Flow, PHX ] = HeatExchanger( Flow, Tout )
eff = 0.95;

H = NetFlow(Flow)*SpecHeat(Flow)*(Flow.T-Tout);
PHX = eff*H/1000;%in kW

Flow.T = Tout;
end

