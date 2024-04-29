function [AutoCorr, TgPosition, FreqOffset] = AutoCorrFunction(RxSignal, WidthWindow, Nfft)
  % Рассчитываем автокорреляционную функцию
  AutoCorr = zeros(1, length(RxSignal)-WidthWindow-Nfft);
  for n = 1:length(RxSignal)-WidthWindow-Nfft
    AutoCorr(n) = sum(RxSignal(n:n+WidthWindow-1).*conj(RxSignal(n+Nfft:n + WidthWindow + Nfft-1)));
    AutoCorr(n) = AutoCorr(n) / sqrt(sum(abs(RxSignal(n:n+WidthWindow-1)).^2)*sum(abs(RxSignal(n+Nfft:n+WidthWindow+Nfft-1)).^2));
  end

  % Поиск позиции защитного интервала с помощью порогового значения
  Threshold = 0.77;
  amp_AutoCorr = abs(AutoCorr);
  Threshold_idxs = find(amp_AutoCorr>Threshold);
  Threshold_idxs = Threshold_idxs(Threshold_idxs>WidthWindow);

  diffs = diff(Threshold_idxs);
  mask = [true, abs(diffs) ~= 1];
  result = find(mask);
  
  try 
    TgPosition=floor((Threshold_idxs(result(1))+Threshold_idxs(result(2)-1))/2);
  catch 
    warning('Не удалось определить позицию защитного интервала. Использовано стандартное значение.');
    TgPosition = 65;
  end

  % Оценка частотной рассинхронизации
  FreqOffset = -angle(AutoCorr(TgPosition))/(2*pi); 
end