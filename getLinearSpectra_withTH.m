function [S,H] = getLinearSpectra_withTH(T, z, ag, Dtg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%The script takes the following input :
%
%T   : Spectral acceleration period for the computation  
%z   : Damping ratio (z = 0.005 for 5% damping)  
%ag  : The acceleration time history of input ground-motion  
%Dtg : The discreet time interval at which the input acceleration is recorded  
%
%The script returns the following:
%
%S  
%|-S.d : Spectral displacement  
%|-S.v : pseudo spectral velocity [S.d * ( 2 * pi / T )]  
%|-S.a : pseudo spectral acceleration [S.d * ( 2 * pi / T ) ^ 2]  
%
%H  
%|-H.d : Time history of the displacement response of the oscillator  
%|-H.v : Time history of the velocity response of the oscillator  
%|-H.a : Time history of the acceleration response of the oscillator  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Interpolate p=-ag*m (linearly)
% ------------------------------
minT = min(T);
Dt = min([minT/10 Dtg 0.55*minT]);
tg = [ 0:(length(ag)-1) ]' * Dtg;
t = [ 0:Dt:tg(end) ]';
m=1;
p = interp1( tg, -ag*m, t );
Gamma = 1/2;
Beta = 1/6; % for linear acceleration
% Beta = 1/4; % for avg acceleration

[sd,d,v,a] = linearNewmark_TH(T,z,p,Gamma,Beta,Dt);
S.d = sd';
S.v = S.d .* (2*3.1415./T);
S.a = S.d .* (2*3.1415./T).^2;

H.d = interp1( t, d', tg );
H.v = interp1( t, v', tg );
H.a = interp1( t, a', tg );