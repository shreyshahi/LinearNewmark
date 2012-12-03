function [S,H] = getLinearSpectra_withTH(T, z, ag, Dtg)

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