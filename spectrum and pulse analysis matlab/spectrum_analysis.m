run('constants.m');                             % Initialize constants (see constants.m for defined constants)
load('power_scale.mat');                        % Load matlab variable power_scale (multiplicative factor to account for detector response)

filename = "20180815 NIR Spectrum at AC.txt";   % Specify filename of spectrum text file
% filename = "2319 spectrum 830gmm.txt";

A = importdata(filename);                       % Import data structure
wavelength = A.data(:,1).*10^-3;                % Collect wavelength in microns from data structure
spectrumRaw = A.data(:,2);                      % Collect raw spectrum in A.U. from data structure
spectrumScaled = spectrumRaw .* power_scale;    % Scale raw spectrum to account for detector response

wavLimits = [wavelength(1), wavelength(end)];   % Find wavelength limits collected by spectrometer
angFreqLimits = 2*pi*c./(wavLimits * 1e-6);     % Convert to angular frequency (Hz)
angFreqCenter = sum(angFreqLimits)/2;           % Calculate center angular frequency (Hz)

%%%% Expand frequency grid for higher time resolution in FFT %%%%
%%%% Pad spectrum to account for expanded angular frequency grid %%%%
angFreq = 2*pi*c./(wavelength * 1e-6);                                      % Angular frequencies calculated from linear wavelength grid
angFreqLinear = linspace(angFreq(1),angFreq(end),length(wavelength));       % Create linear-spaced angular frequency grid
deltaAngFreq = abs(angFreqLinear(1) - angFreqLinear(2));                    % Calculate delta angular frequency

N = 2^14;                                                                   % Number of points to use in FFT
spectrumAngFreq = spline(angFreq, spectrumScaled, angFreqLinear);           % Spline to calculate spectrum for linear angular frequency grid
numPad = (N - length(spectrumAngFreq))/2;                                   % Calculate number of zeros to pad on front and back of spectrum
spectrumAngFreqExpanded = padarray(spectrumAngFreq,[0,numPad],0,'both');    % Zero-padded spectrum
angFreqLinearHigh = fliplr(angFreqLinear(1) + deltaAngFreq.*(1:numPad));    % angular frequencies to pad on high frequency side
angFreqLinearLow = angFreqLinear(end) - deltaAngFreq.*(1:numPad);           % angular frequencies to pad on low frequency side
angFreqExpanded = [angFreqLinearHigh, angFreqLinear, angFreqLinearLow];     % Concatonate angular frequency arrays

%%%% Spectral Filter to reduce noise or filter spectrum %%%%
spectralFilterWav = [0.65 1.05];
% spectralFilterWav = [0.69 0.7];
spectralFilterAngFreq = 2*pi*c./(spectralFilterWav * 1e-6);
spectrumAngFreqExpanded(angFreqExpanded > spectralFilterAngFreq(1)) = 0;
spectrumAngFreqExpanded(angFreqExpanded < spectralFilterAngFreq(2)) = 0;

figure(1); plot(angFreqExpanded, spectrumAngFreqExpanded)

%%%% Setup FFT for calculate transform-limited pulse duration %%%%
deltaTime = 2*pi/(N*deltaAngFreq);      % Calculate delta time (s)
time = (-N/2:N/2-1)*deltaTime;          % Create time grid

ampAngFreq = sqrt(spectrumAngFreqExpanded);                                             % Cacluate amplitude of E-field as a function of angular freuqency;
ampTime = fftshift(ifft(ifftshift(ampAngFreq .* exp(-1i * angFreqCenter * time))));     % IFFT to calcuate E-field in time domain
%%%% Note about shift in frequency grid: exp(1i * angFreqCenter * time) needed to plot. No need to shift again because already applied once
intensity = ampTime.*conj(ampTime);

figure(2); plot(time,ampTime.*conj(ampTime), 'k');      % Plot intensity in time domain to check pulse
T_duration = pulse_duration(time, intensity, 2);        % Calculate pulse duration. See pulse_duration.m for details.
fprintf('FWHM Pulse Duration: %0.2e s\n', T_duration)

%%%% Calculate autocorrelation trace %%%%
AC = xcorr(intensity,intensity);                        % Calculate autocorrelation of pulse
numPoints = length(AC);                                 % Check number of points in autocorrelation to create new time grid

timeCorr = (time(2)-time(1)) * (-(numPoints-1)/2:(numPoints-1)/2);  % create new time grid with numPoints of points

figure(3); plot(timeCorr,AC,'r');
AC_duration = pulse_duration(timeCorr, AC, 2);
fprintf('AC FWHM: %0.2e s\n', AC_duration)

%%%% Calculate autocorrelation by fourier transformation of spectral power
%%%% density squared
ampSpectrumFromFFT = ifftshift(fft(fftshift(ampTime))) .* exp(1i * angFreqCenter * time);     % IFFT to calcuate E-field in time domain
%%%% Note about shift in frequency grid: exp(1i * angFreqCenter * time)
%%%% applied becasue wanted to plot spectrum
figure(4); plot(angFreqExpanded, spectrumAngFreqExpanded, angFreqExpanded, ampSpectrumFromFFT.*conj(ampSpectrumFromFFT),'r--')

intensityFFT = ifftshift(fft(fftshift(intensity))) ;        % FFT of intensity in time domain to calculate AC
AC_intensity = fftshift(ifft(ifftshift(intensityFFT.^2)));  % AC calculated from IFFT( FFT(intensity).^2 )

[~,maxID] = max(AC_intensity);
figure(5); plot(time - time(maxID), AC_intensity, timeCorr, AC);
% figure(5); plot(time, AC_intensity, timeCorr, AC);

%%%% Apply dispersion to try to match uncompressed autocorrelation
centerWavelengthDispersion = 0.761;
centerAngFreqDispersion = 2*pi*c./(centerWavelengthDispersion .* 1e-6);

GDD = 1800 * 1e-30;
TOD = -4000 * 1e-45;

ampAngFreqStreteched = ampAngFreq .* exp(1i * (...
    GDD/2 * (angFreqExpanded - centerAngFreqDispersion).^2 + ...
    TOD/6 * (angFreqExpanded - centerAngFreqDispersion).^3));
ampTimeStretched = fftshift(ifft(ifftshift(ampAngFreqStreteched .* exp(-1i * angFreqCenter * time))));
intensityStretched = ampTimeStretched.*conj(ampTimeStretched);

%%%% IFFT( FFT(intensity).^2 ) doesn't seem to work well. Perhaps run out
%%%% of time grid on the right hand side.
% intensityStretchedFFT = ifftshift(fft(fftshift(intensityStretched))) ;        % FFT of intensity in time domain to calculate AC
% AC_intensityStretched = fftshift(ifft(ifftshift(intensityStretchedFFT.^2))) .* exp(1i * angFreqCenter * time);  % AC calculated from IFFT( FFT(intensity).^2 )
% 
% [~,maxIDStretched] = max(AC_intensityStretched);
% figure(6); plot(time - time(maxIDStretched), AC_intensityStretched);

%%%% Calculate autocorrelation trace %%%%
ACStretched = xcorr(intensityStretched,intensityStretched);                        % Calculate autocorrelation of pulse
numPoints = length(ACStretched);                                 % Check number of points in autocorrelation to create new time grid

timeCorrStretched = (time(2)-time(1)) * (-(numPoints-1)/2:(numPoints-1)/2);  % create new time grid with numPoints of points

figure(7); plot(timeCorrStretched,ACStretched/max(ACStretched),'r');
AC_durationStretched = pulse_duration(timeCorrStretched, ACStretched, 2);
fprintf('AC FWHM: %0.2e s\n', AC_durationStretched)