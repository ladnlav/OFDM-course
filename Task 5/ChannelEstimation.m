pskModulator = comm.PSKModulator;
pskModulator.ModulationOrder = 4;
N = 256;
K = 4;
insig = randi([0,pskModulator.ModulationOrder-1],N,1);
constpoints = pskModulator(insig);
channelInput = ifft(constpoints, K*N);

ricianchan = comm.RicianChannel( ...
    'SampleRate',1e6, ...
    'PathDelays',[0.0 200]*1e-6, ...
    'AveragePathGains',[1 20], ...
    'KFactor',2.8, ...
    'DirectPathDopplerShift',0.0, ...
    'DirectPathInitialPhase',0.0, ...
    'MaximumDopplerShift',0, ...
    'DopplerSpectrum',doppler('Bell', 8), ...
    'RandomStream','mt19937ar with seed', ...
    'Seed',73, ...
    'PathGainsOutputPort',true);

[RicianChanOut1,RicianPathGains1] = ricianchan(channelInput);

H = RicianChanOut1 .* conj(channelInput);
plot(10*log10(abs(ifft(H))))