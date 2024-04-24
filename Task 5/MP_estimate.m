% MP оценка
function [H_MP,h_impulse_est] = MP_estimate(Y, pilot_indices, Nfft, N_carrier,dominant_taps)
Np = length(pilot_indices);
comb = pilot_indices(3)-pilot_indices(2);

h_impulse_est = zeros(Nfft,1);
residue = Y(pilot_indices,1);
F = dftmtx(Nfft); 
F = F(:,1:ceil(Nfft/comb));
P = zeros(Np, Nfft); 
for i=1:Np
    P(i,pilot_indices(i)) = 1;
end
sensing_matrix = P*F;
kp=zeros(1,dominant_taps);
x = zeros(1, dominant_taps);

for i1 = 1:dominant_taps
    projection = 0;
        for i2 = 1:Np 
            if ismember(i2,kp)
                projection(i2)=-100;
            else
                a=sensing_matrix(:,i2);
                projection(i2) = abs(a'*residue)^2/(norm(a)^2);
            end
        end
    [~,kp(i1)] = max(projection);
    
    % residual vector
    proj_kp = sensing_matrix(:,kp(i1))*sensing_matrix(:,kp(i1))'/(norm(sensing_matrix(:,kp(i1)))^2);
    x(i1) = sensing_matrix(:,kp(i1))'*residue/(norm(sensing_matrix(:,kp(i1)))^2);
    residue = residue - proj_kp*residue;
end

kp = sort(kp) - min(kp)+1;
h_impulse_est = zeros(Nfft,1);
for i1 = 1:length(kp)
    h_impulse_est(kp(i1)) = x(i1);
end


h_impulse_est = h_impulse_est/1.7889;
H_MP = fft(h_impulse_est).';
end