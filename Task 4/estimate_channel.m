function [H_est,Hest_at_pilots] = estimate_channel(rx_signal,allCarriers,pilotCarriers,pilotValues)
tx_pilotValues = pilotValues;

rx_pilotValues = rx_signal(pilotCarriers,:);

Hest_at_pilots = mean(rx_pilotValues ./ tx_pilotValues,2);

H_est= interp1(pilotCarriers, Hest_at_pilots, allCarriers, 'spline');

end

