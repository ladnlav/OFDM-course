% MP оценка
function [H_OMP,h_impulse_est,index] = OMP_estimate(Y, pilot_indices, Nfft, N_carrier,dominant_taps)

Np = length(pilot_indices);
comb = pilot_indices(3)-pilot_indices(2);

h_impulse_est = zeros(Nfft,1);
measurement_vec = Y(pilot_indices,1);
F = dftmtx(Nfft);
F = F(:,1:ceil(Nfft/comb));
P = zeros(Np, Nfft);
for i=1:Np
    P(i,pilot_indices(i)) = 1;
end
sensing_matrix = P*F;
kp=zeros(1,dominant_taps);

index = [];
[~,index(1)] = max(abs(sensing_matrix'*measurement_vec));
A = sensing_matrix(:,index);
x = pinv(A)*measurement_vec; %least squares solution
residue = [];
residue(:,1) = measurement_vec - A*x;

for i1= 2:dominant_taps
    [~,index(i1)] = max(abs(sensing_matrix'*residue(:,i1-1)));

    A = [A sensing_matrix(:,index(i1))];
    x = pinv(A)*measurement_vec;
    residue(:,i1) = measurement_vec - A*x;
end

len = length(x);
x = [x;zeros(Nfft-len,1)];
est_fade_chan = zeros(Nfft,1);

index = sort(index);
index = index - min(index)+1;
for i1 = 1:length(index)
    est_fade_chan(index(i1)) = x(i1);
end
h_impulse_est = est_fade_chan.'/1.7889; % now a row vector

H_OMP = fft(h_impulse_est);
end