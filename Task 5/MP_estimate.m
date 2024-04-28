% MP оценка
function [H_MP,h_impulse_est] = MP_estimate(Y, sensing_matrix, Nfft,dominant_taps)
Np = size(sensing_matrix,1);
residue = Y;
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

% kp = sort(kp) - min(kp)+1;
h_impulse_est = zeros(Nfft,1);
for i1 = 1:length(kp)
    h_impulse_est(kp(i1)) = x(i1);
end

h_impulse_est = h_impulse_est;
H_MP = fft(h_impulse_est).';
end