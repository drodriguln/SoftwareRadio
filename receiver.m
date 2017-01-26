clc;
clear all;
close all;

%%-- File Read --%%

[st,Fs,nbits] =  wavread('modulate.wav');

N = length(st);
f = (Fs/N)*[0:1:N-1];
df = Fs/N;
Ts = 1/Fs;
t = Ts*[0:1:N-1];

%%-- Bandpass Demodulation --%%

[Hzn_BP1, Hzd_BP1] = butter(2, [5000 13000]/(Fs/2),'bandpass');
demod1 = filter(Hzn_BP1,Hzd_BP1,st);
 
[Hzn_BP2, Hzd_BP2] = butter(2, [15000 23000]/(Fs/2),'bandpass');
demod2 = filter(Hzn_BP2,Hzd_BP2,st);

demod1_fft = abs(fft(demod1));
demod2_fft = abs(fft(demod2));

figure(1)
plot(f, demod1_fft)
xlim([0 Fs/2])
grid on

figure(2)
plot(f, demod2_fft)
xlim([0 Fs/2])
grid on

%%-- Local Oscillator --%%

fc1 = 9000;
fc2 = 19000;

% Synchronization of the demod1's Local Oscillator
k = 0;
threshold = 0.5*max(st);
for i=2:Fs % about 1 sec of s(t)
   if (demod1(i-1) < demod1(i) & demod1(i) >= ...
           demod1(i+1) & demod1(i) > threshold)
         k = k + 1;
         peaks(k) = i;
   end;
end;
phi1  = 2*pi*fc1*t(peaks(1));
demod1_LO = cos(2*pi*fc1*t - phi1);
% -------------------------------------------------------------

% Synchronization of the demod2's Local Oscillator
k = 0;
threshold = 0.5*max(st);
for i=2:Fs % about 1 sec of s(t)
   if (demod2(i-1) < demod2(i) & demod2(i) >= ...
           demod2(i+1) & demod2(i) > threshold)
         k = k + 1;
         peaks(k) = i;
   end;
end;
phi2  = 2*pi*fc2*t(peaks(1));
demod2_LO = cos(2*pi*fc2*t - phi2);
% -------------------------------------------------------------

st1 = demod1.*demod1_LO';
st2 = demod2.*demod2_LO';

%%-- Low-Pass Filter --%%

f0 = 4000;
wn = f0/(Fs/2);

[Hzn_LP, Hzd_LP] = butter(10, wn);

st1_LP = filter(Hzn_LP,Hzd_LP,st1);
st2_LP = filter(Hzn_LP,Hzd_LP,st2);

figure(3)
plot(t, st1_LP)
grid on

figure(4)
plot(t, st2_LP)
grid on

%%-- DC Offset Removal --%%

Ac = 0.25;

mt1 = st1_LP - Ac/2;
mt2 = st2_LP - Ac/2;

figure(5)
plot(t, mt1)
grid on

figure(6)
plot(t, mt2)
grid on

figure(7)
plot(f, abs(fft(mt1)))
xlim([0 Fs/2])
grid on

figure(8)
plot(f, abs(fft(mt2)))
xlim([0 Fs/2])
grid on

%%-- File Write --%%

wavwrite(mt1,Fs,nbits,'mt1.WAV')
wavwrite(mt2,Fs,nbits,'mt2.WAV')