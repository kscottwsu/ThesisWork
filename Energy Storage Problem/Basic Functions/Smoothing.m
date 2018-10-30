function [ SProfile ] = Smoothing( Profile, points, value )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

SProfile = 0*Profile;
LenProfile = length(Profile);

for i = 1:LenProfile
    if i > points && i+points < LenProfile  %% not near bounds
        SProfile(i) = sum(Profile((i-points):(i+points)))/(2*points+1);
    elseif i <= points && i+points < LenProfile  %% left bound
        Lbound = LenProfile + (i -1) - points;
        SProfile(i) = (sum(Profile(Lbound:end)) + sum(Profile(1:(i+points))))/(2*points+1);
    elseif i > points && i+points >= LenProfile %% right bound
        Rbound = i+points-LenProfile;
        SProfile(i) = (sum(Profile(1:Rbound)) + sum(Profile((i-points):LenProfile)))/(2*points+1);
    end
end

SProfile = SProfile*value;
end

