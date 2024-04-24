% Рассчитать CCDF
function [PAPR_ccdf, CCDF] = calculateCCDF(PAPR_values)
    
    [CCDF,PAPR_ccdf] = ecdf(PAPR_values);
    CCDF = 1 - CCDF;
end