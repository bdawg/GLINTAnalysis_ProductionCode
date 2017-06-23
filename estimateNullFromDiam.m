%B=6.59; 
B=5.65; %Assumes 5.10 mm pupil at MEMS, 7.92m primary pupil
B=5.55; %Assumes 5.10 mm pupil at MEMS, 7.79m primary pupil (ie after scexao mask)

lambda = 1.55e-6;
ldcoeff = 0.;
theta = 16.3; %mas diameter

%ldcoeff = 0.35;
theta = 31.5; %mas diameter

% lambda=2.16e-6;
% B=3.2;
% theta=20.95;


thetaRad = theta/1000/60/60/360*2*pi;
% N = ( (pi*B*thetaRad)/(4*lambda) )^2 * (1 - ((7*ldcoeff)/15)) * ...
%     (1 - ldcoeff/3)
N = ( (pi*B*thetaRad)/(4*lambda) )^2 * (1 - ((7*ldcoeff)/15)) * ...
    (1 - ldcoeff/3)^(-1) %The ^-1 is in Absil2011 but not Hanot!