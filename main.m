numSymbs = 2^16;
M = 4; % QPSK

Rsym = 28e9; % symbol rate (sym/sec)

%% zs: array of distances z to be simulated
% Example: zs = 42; zs = 40:10:100; zs = [300, 500, 1000];
zs = 42;
plotlen = length(zs);

%% Tx RRC filter properties
rolloff = 0.25;
span = 6; % filter span
sps = 16; % samples per symbol

%% Sampling frequency
fs = Rsym * sps; % Hz
Tsamp = 1 / fs; % s
% t: time vector, s
t = (0 : 1 / fs : numSymbs / Rsym - 1 / fs).';


%% Chromatic dispersion
D = 17; % ps / (nm km)
lambda = 1550; % nm

%% Laser phase noise
linewidthTx = 0; % Hz
linewidthLO = 1e6; % Hz

%% Kerr effect / SSFS parameters
gamma = 1.2; % watt^-1 / km
alpha = 0.2; % dB/km
dz = 2; % Step size, km

%% Polarization state rotation parameters
rot_omega = 1e3; % rad/s
rot_phi = 2; % rad

%% Launch power, per wavelength channel
power_dBm = 0;
power = 10 .^ (power_dBm / 10) * 1e-3; % watts

%% WDM properties
wavelength_channels = 3;
dw = 2 * pi * 50e9; % channel spacing (rad/s)

%% Shot noise
hc = 6.62607015e-34 * 299792458; % J m
Eperphoton = hc / (lambda * 1e-9); % J


%% Stores result to be plotted
ber = zeros(plotlen, 1);
if plotlen > 1
  fig = figure;
end

%% sps and Tsamp change at Tx/Rx, save these for later.
spsOrig = sps;
TsampOrig = Tsamp;

%% Generate random data for both polarizations
data_x = randi([0, M - 1], numSymbs, wavelength_channels, 'uint8');
data_y = randi([0, M - 1], numSymbs, wavelength_channels, 'uint8');

%% DE-QPSK modulation
modData_x = deqpskmod(data_x);
modData_y = deqpskmod(data_y);

%% Construct waveforms for each channel separately
A_x_wdm = zeros(numSymbs * sps, wavelength_channels);
A_y_wdm = zeros(numSymbs * sps, wavelength_channels);
carriers = zeros(numSymbs * sps, wavelength_channels);

for w = 1 : wavelength_channels
  %% Compute frequency offsets:
  %                     ___     ___     ___     ___     ___
  %       Spectrum     |   |   |   |   |   |   |   |   |   |
  %                    |   |   |   |   |   |   |   |   |   |
  %                ____|   |___|   |___|   |___|   |___|   |____  --> freq
  % channel #            5       3       1       2       4
  % ang freq offset    -2dw     -dw      0      +dw    +2dw

  if mod(w, 2) == 0
    ndw = w / 2 * dw;
  else
    ndw = (1-w) / 2 * dw;
  end
  carriers(:, w) = exp(1j * ndw * t);
  A_x_wdm(:, w) = txFilter(modData_x(:, w), rolloff, span, sps);
  A_y_wdm(:, w) = txFilter(modData_y(:, w), rolloff, span, sps);
end

%% Sum the WDM waveforms with their frequency offsets
A_x = sum(A_x_wdm .* carriers, 2);
A_y = sum(A_y_wdm .* carriers, 2);

%% Clear variables no longer needed to reduce memory usage
clear modData_x modData_y A_x_wdm A_y_wdm;

%% Set launch power. Divide by 2 because half power for each polarization.
A_x = sqrt(power / 2) * A_x;
A_y = sqrt(power / 2) * A_y;

%% Rotate polarization states
A_x = A_x .* cos(rot_omega * t) + ...
      A_y .* sin(rot_omega * t) * exp(-1j * rot_phi);
A_y = A_x .* -sin(rot_omega * t) * exp(1j * rot_phi) + ...
      A_y .* cos(rot_omega * t);

%% Now loop through each z
for i = 1 : plotlen
  z = zs(i);

  sps = spsOrig;
  Tsamp = TsampOrig;

  %% Split-step Fourier
  [Al_x, Al_y] = ssfs(A_x, A_y, D, lambda, z, dz, Tsamp, gamma, alpha);

  %% Phase noise
  Al_x = phaseNoise(Al_x, linewidthTx, linewidthLO, Tsamp);
  Al_y = phaseNoise(Al_y, linewidthTx, linewidthLO, Tsamp);

  %% Here, only receive the central channel 1.
  % For channel n: Al_x .* conj(carriers(:, n)); etc.
  r_x = rxFilter(Al_x, rolloff, span, sps);
  r_y = rxFilter(Al_y, rolloff, span, sps);
  % Rx filter performs downsampling as well, keep track of this
  sps = 2;
  Tsamp = Tsamp * spsOrig / sps;

  %% Rx shot noise
  photonpersym = mean(abs(r_x) .^ 2) / Rsym / Eperphoton;
  snr = photonpersym;
  r_x = awgn(r_x, snr, 'measured', 'linear');
  r_y = awgn(r_y, snr, 'measured', 'linear');

  %% -- Begin DSP channel equalization --
  %% Chromatic dispersion compensation
  r_x = CDCompensation(r_x, D, lambda, z, Tsamp);
  r_y = CDCompensation(r_y, D, lambda, z, Tsamp);
  r_x = r_x(2:2:end);
  r_y = r_y(2:2:end);

  %% Adaptive filter
  [r_x, r_y] = pdm_adaptiveCMA(r_x, r_y);

  %% Phase noise correction
  r_x = phaseNoiseCorr(r_x, M, 0, 40).';
  r_y = phaseNoiseCorr(r_y, M, 0, 40).';

  %% Demodulate DE-QPSK
  demod_x = deqpskdemod(r_x);
  demod_y = deqpskdemod(r_y);

  %% Calculate and store BER
  [~, ber(i)] = biterr([data_x(:, 1); data_y(:, 1)], [demod_x; demod_y]);

  q = 20 * log10(erfcinv(2*ber)*sqrt(2));
  if i > 1
    figure(fig);
    plot(zs, q);
  end
end

ber
q
