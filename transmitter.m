clc;
clear;
close all;

%%-- File Read --%%

[station1,Fs,nbits] =  wavread('Soundtrack-Theme_12.wav');
[station2,Fs,nbits] =  wavread('Bloody_Tears.wav');

%%-- File Vector Length Fix --%%

N2 = length(station2);
N1 = length(station1);
station1_fix = station1(1:N2);

%%-- Low-Pass Filter --%%

f0 = 4000;
df = Fs/N2;
Ts = 1/Fs;
wn = f0/(Fs/2);
t = Ts*[0:1:N2-1];
f = (Fs/N2)*[0:1:N2-1];

[Hzn_LP, Hzd_LP] = butter(10, wn);

mt1 = filter(Hzn_LP,Hzd_LP,station1_fix);
mt2 = filter(Hzn_LP,Hzd_LP,station2);

%%-- AM Modulation --%%

Ac = 0.25;
fc1 = 9000;
fc2 = 19000;

st1 = Ac*[1 + mt1].*cos(2*pi*fc1*t');
st2 = Ac*[1 + mt2].*cos(2*pi*fc2*t');

%%-- AM Signal Combination --%%

st = st1 + st2;
st_fft = fft(st);

figure(1)
plot(mt1)
grid on;

figure(2)
plot(st1)
grid on;

figure(3)
plot(f, abs(st_fft))
xlim([0 (Fs/2)])
grid on;

%%-- File Write --%%

wavwrite(st,Fs,nbits,'modulate.WAV')
