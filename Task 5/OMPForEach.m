 function [h_t_s, h_f_s] = OMPForEach(Y, Nfft, pilot_ind, max_it)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

h_f = Y;
h_f_ = Y(pilot_ind,1);
Np = length(pilot_ind);

F = dftmtx(Nfft);
h_t_s = zeros(Nfft, max_it);
h_f_s = zeros(Np, max_it);
% corr_s = zeros(1, max_it);
indexes = zeros(1, max_it);
res_i = h_f_;
PF_ = F(pilot_ind, :);

for it = 1:max_it
    it;
%     res_prev = res_i;

    [~, indexes(it)] = max(abs(res_i' * PF_));
    h_t = pinv(PF_(:, indexes(1:it))) * h_f_;
    res_i = h_f_ - PF_(:, indexes(1:it)) * h_t;
    h_t_s(indexes(1:it),it) = h_t;
    h_f_s(:,it) = PF_(:, indexes(1:it))*h_t;

    PF_(:,indexes(it)) = 0;
%     corr_s(it) = abs(h_f' * h_f_s(:,it))/norm(h_f)/norm(h_f_s(:,it));
    
%     if (norm(res_i) < sqrt(noise_pw)) || (norm(res_i - res_prev) / norm(res_prev) < 1e-2)
%         break;
%     end
end

end