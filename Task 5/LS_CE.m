function [H_LS] = LS_CE(Y,Xp,pilot_loc,Nfft)
% LS channel estimation function
% Inputs:
% Y = Frequency-domain received signal
% Xp = Pilot signal
% pilot_loc = Pilot location
% N = FFT size
% Nps = Pilot spacing
% int_opt = ’linear’ or ’spline’
% output:
% H_LS = LS Channel estimate


Np=length(pilot_loc);   
% 
% F = dftmtx(Nfft); 
% P = zeros(Np, Nfft); 
% for i=1:Np
%     P(i,pilot_loc(i)) = 1;
% end
% sensing_matrix = P*F;
% 
% h_t_ls = pinv(sensing_matrix)*Y(:);
% H_LS = (h_t_ls).';

% h_t_ls = pinv(sensing_matrix)*Y(pilot_loc,1);
% H_LS = fft(h_t_ls).';


k=1:Np;
LS_est(k) = Y(pilot_loc(k))./Xp(k); % LS channel estimation

% Linear/Spline interpolation
H_LS = interpolate(LS_est,pilot_loc,Nfft,'spline');


end