% MP оценка
function [H_OMP,h_impulse_est,index] = OMP_estimate(Y, sensing_matrix, Nfft,dominant_taps, SNR_dB)

noise_pw = 10^(SNR_dB/10);
measurement_vec = Y;
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

    if ((norm(residue(:,i1)-residue(:,i1-1))/norm(residue(:,i1-1)) < 1e-2))
            break;
    end
end

len = length(x);
x = [x;zeros(Nfft-len,1)];
est_fade_chan = zeros(Nfft,1);

% index = sort(index);
% index = index - min(index)+1;
for i1 = 1:length(index)
    est_fade_chan(index(i1)) = x(i1);
end
h_impulse_est = est_fade_chan.'; % now a row vector

H_OMP = fft(h_impulse_est);
end