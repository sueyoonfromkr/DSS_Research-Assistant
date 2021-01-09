
% name: Sue Yoon
% date: Aug 26 2020

clear;
eeglab;

data = ALLEEG.data;


%% create a complex Morlet wavelet (delta)

srate = 1200;
time  = -1.5:1/srate:1.5;
frex = 2.5;

sine_wave = exp( 1i*2*pi*frex*time );
fwhm = 0.8; % width of the Gaussian in seconds
gaus_win = exp( (-4 * log(2)*time.^2) / fwhm^2 );

wavelet = sine_wave .* gaus_win;
half_wavN = (length(time)-1)/2;



% plot the wavelet in the time domain
figure(2)
plot(time, real(wavelet));

% plot the wavelet in the frequency domain
pnts = length(time);

mwX = abs(fft( wavelet )/pnts);
hz  = linspace(0,srate,pnts);

figure(3)
plot(hz,mwX,'k','linew',2)
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Complex Morlet wavelet in the frequency domain')






%% convolution

% FFT parameters
nWave = length(time);
nData = size(data,2);
nConv = nWave + nData - 1;

% FFT of wavelet (check nfft)
waveletX = fft(wavelet,nConv);
waveletX = waveletX ./ max(waveletX);

figure(4)
plot(abs(waveletX));





%% phase synchronization

% initialize output time-frequency data
timevec = EEG.times;
phase_data = zeros(2, length(timevec)); 
real_data  = zeros(2, length(timevec));

    %%% for loop should be done here!!! 
chan1idx = 3;
chan2idx = 14;

% analytic signal of channel 1
dataX = fft(data(chan1idx, :), nConv);
as = ifft(waveletX.*dataX,nConv);
as = as(half_wavN+1:end-half_wavN);

% collect real and phase data
phase_data(1,:) = angle(as); % extract phase angles
real_data(1,:)  = real(as);  % extract the real part (projection onto real axis)
plot(timevec, phase_data);
plot(timevec, real_data);

% analytic signal of channel 8
dataX = fft(data(chan2idx, : ), nConv);
as = ifft(waveletX.*dataX,nConv);
as = as(half_wavN+1:end-half_wavN);

% collect real and phase data
phase_data(2,:) = angle(as);
real_data(2,:)  = real(as);

figure(5)
subplot(211)
plot(timevec, phase_data);

subplot(212)
plot(timevec, real_data);


%%  quantify phase synchronization between the two channels

% phase angle differences
% on one line:
phase_synchronization = abs(mean(exp(1i * diff(phase_data))));

disp([ 'Synchronization between ' num2str(chan1idx) ' and ' num2str(chan2idx) ' is ' num2str(phase_synchronization) '!' ])

% on one line:
phase_synchronization = abs(mean(exp(1i * diff(phase_data))));


%% for loop per channels 

timevec = EEG.times;
phase_data = zeros(2, length(timevec)); 
real_data  = zeros(2, length(timevec));


band = zeros(5, 91)
band_idx = 1

frex_list = [2.5 6 10.5 21.5 40];
fwhm_list = [0.8 0.5 0.3 0.13 0.27];


for s = 1:5
    srate = 1200;
    time = -1.5: 1/srate : 1.5;
    frex = frex_list(s);
    fwhm = fwhm_list(s);
    
    sine_wave = exp( 1i*2*pi*frex*time );
    gaus_win = exp( (-4 * log(2)*time.^2) / fwhm^2 );

    wavelet = sine_wave .* gaus_win;
    half_wavN = (length(time)-1)/2;
    
    nWave = length(time);
    nData = size(data,2);
    nConv = nWave + nData - 1;
    
    waveletX = fft(wavelet,nConv);
    waveletX = waveletX ./ max(waveletX);
    
    band_idx = 1
    
    for a = 1:13
        chan1idx = a;
        for j = a+1 : 14
            chan2idx = j;
        
            dataX = fft(data(chan1idx, :), nConv);
            as = ifft(waveletX.*dataX,nConv);
            as = as(half_wavN+1:end-half_wavN);
        
            phase_data(1,:) = angle(as); % extract phase angles
            real_data(1,:)  = real(as);  % extract the real part (projection onto real axis)

        
            dataX = fft(data(chan2idx, : ), nConv);
            as = ifft(waveletX.*dataX,nConv);
            as = as(half_wavN+1:end-half_wavN);
        
            phase_data(2,:) = angle(as);
            real_data(2,:)  = real(as);
        
            band(s, band_idx) = abs(mean(exp(1i * diff(phase_data))));
        
            band_idx = band_idx + 1;


        
        end
    end
end


band


subject_num = {EEG.setname(1); EEG.setname(1); EEG.setname(1); EEG.setname(1) ; EEG.setname(1)};
band_name = {'delta'; 'theta'; 'alpha'; 'beta'; 'gamma'};
ses = {'hri'; 'hri'; 'hri'; 'hri'; 'hri'};
band = [subject_num, band_name, ses, num2cell(band)];

file_name = strcat('ps_', EEG.setname(1), 'hri.csv');

writecell(band, file_name)


