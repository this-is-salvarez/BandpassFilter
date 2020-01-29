clear all; close all;

%Parameters for Filter
L = 550;
M = L/2;
Z = L+1;
N = 80000; 
Fs = 8000; 
Fh = 320+15;
f_H = Fh/Fs; 
Fl = 320-15;
f_L = Fl/Fs; 
freq = (0:N/2)*Fs/N;
f_plot = (-N/2:(N/2-1))*Fs/N; 
%Create Empty Vectors to be filled in the for loop
h=zeros(Z,1);
hL=zeros(Z,1);
hH=zeros(Z,1);
Hann_window = zeros(Z,1);
Hamming_window = zeros(Z,1);
Blackman_window = zeros(Z,1);

for n = 1:Z
    hH(n) = 2*f_H*sinc(2*f_H*(n-M-1));
    hL(n) = 2*f_L*sinc(2*f_L*(n-M-1));
    h(n) = (hH(n) - hL(n));
    %Best window. Has better attenuation at stopband
    Hann_window(n) = 1/2*(1-cos(2*pi*(n)/(Z-1)));
    %Other windows
    Hamming_window(n) = (0.54 - 0.46*cos((2*pi*n)/Z-1));
    % Hann_window(n)= 1/2*(1-cos(2*pi*(n)/(Z-1)));
    Blackman_window(n)=0.42-0.5*cos(2*pi*n/(Z-1))+0.08*cos(4*pi*n/(Z-1));
end
h_Hann = Hann_window.*h;
h_Hamming = Hamming_window.*h;
h_Blackman = Blackman_window.*h;

% %Plot the Magnitude and Phase for the Lower Low Pass Filter
% HL1 = fft(hL,N);
% HL1_Mag = (abs(HL1(1:N/2+1)));
% HL1_angle = angle(HL1(1:N/2+1));
% figure(1);
% plot(f_plot',fftshift(abs(HL1)));
% xlabel('Frequency (Hz)'); 
% ylabel('Magnitude - |X(f)|');
% title('Magnitude plot Lower Low Pass Fitler');
% figure(2);
% plot(freq,HL1_angle);
% xlabel('Frequency (Hz)'); 
% ylabel('Phase - <X(f) (rads)');
% title('Phase plot Lower Low Pass Fitler');

% %Plot the Magnitude and Phase for the Higher Low Pass Filter
% HL2 = fft(hH,N);
% HL2_Mag = abs(HL2(1:N/2+1));
% HL2_angle = angle(HL2(1:N/2+1));
% figure (3);
% plot(f_plot',fftshift(abs(HL2)));
% xlabel('Frequency (Hz)'); 
% ylabel('Magnitude - |X(f)|');
% title('Magnitude plot Higher Low Pass Fitler');
% figure(4);
% plot(freq,HL2_angle);
% xlabel('Frequency (Hz)'); 
% ylabel('Phase - <X(f) (rads)');
% title('Phase plot Higher Low Pass Fitler');

% %Plot the Magnitude and phase of the Filter Function no filter
% H = fft(h,N);
% H_Mag = abs(H(1:N/2+1));
% H_angle = angle(H(1:N/2+1));
% figure(5);
% plot(f_plot',fftshift(abs(H)));
% xlabel('Frequency (Hz)'); 
% ylabel('Magnitude - |X(f)|');
% title('Magnitude plot Bandpass Fitler');
% figure(6);
% plot(freq,H_angle);
% xlabel('Frequency (Hz)'); 
% ylabel('Phase - <X(f) (rads)');
% title('Phase plot Bandpaass Filter');

% %Plot the Magnitude Window Filter
H_Hann = fft(h_Hann,N);
H_Hamming =  fft(h_Hamming,N);
H_Blackman =  fft(h_Blackman,N);
H_Hann_Mag = abs(H_Hann(1:N/2+1));
H_Hann_angle = angle(H_Hann(1:N/2+1));
figure(7);
plot(f_plot',fftshift(abs(H_Hann)));
xlabel('Frequency (Hz)'); 
ylabel('Magnitude - |X(f)|');
title('Magnitude plot Bandpass Filter with Hann Window');
figure(8);
plot(f_plot',fftshift(10*log10(abs(H_Hann))));
% hold on 
% plot(f_plot',fftshift(10*log10(abs(H_Hamming))), 'r');
% hold on
% plot(f_plot',fftshift(10*log10(abs(H_Blackman))), 'b');
xlabel('Frequency (Hz)'); 
ylabel('Magnitude -  dB');
xlim([0 640]);
title('Magnitude plot Bandpass Filter with Hann Window ');

figure(9);
plot(freq,H_Hann_angle);
xlabel('Frequency (Hz)'); 
ylabel('Phase - <X(f) (rads)');
title('Phase plot Bandpass Filter with Hann Window');
%Save the filter
 save('BPFilter.mat');