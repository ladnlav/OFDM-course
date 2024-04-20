function y_STO=add_STO(y, nSTO)
    % add STO (symbol time offset)
    % y : Received signal
    % nSTO : Number of samples corresponding to STO
    if nSTO>=0
        y_STO=[y(nSTO+1:end); zeros(1,nSTO)']; % advance
    else 
        y_STO=[zeros(1,-nSTO)'; y(1:end+nSTO)]; % delay
    end
end