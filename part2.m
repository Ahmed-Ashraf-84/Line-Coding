close all;
fc = 10e9;
bit_rate=0.5*10e9;
Tb=1/bit_rate;
bits = randi([0, 1], 1, 100);
n = 100*fc;
T = length(bits)/bit_rate;
N = n * (length(bits)/bit_rate);
dt = T/N;
t = 0:dt:(T - dt);
sig = zeros(1, length(t));
for i = 1:length(bits)
    if bits(i) == 1
        sig(((i - 1) * n)/bit_rate + 1:(i * n)/bit_rate) = 1;
    else
        sig(((i - 1) * n)/bit_rate + 1:(i * n)/bit_rate) = -1;
    end
end

# Spectral Domain
df = 1/T;
fs = 1/ dt;
N = length(t);
if (rem(N,2)==0)
f=- (0.5 * fs) : df : (0.5*fs - df) ;
else
f=- (0.5*fs-0.5*df) : df : (0.5*fs -0.5*df) ;
end

signalFreq = abs(fftshift(fft(sig))/N);
power_signal = signalFreq.^2;
figure(2);
plot(f, power_signal);
axis([0, fs/100, 0, 1.5 * max(power_signal)]);
xlabel('Time');
ylabel('Power');
title('spectral domain');
grid on;

#BPSK Signal

mod_signal = sqrt(2 / Tb) * sig.* cos(2 * pi * fc * t);

figure(3);
plot(t, mod_signal), axis([0, T/10, -1.5*max(mod_signal), 1.5*max(mod_signal)]);
xlabel('Time');
ylabel('Amplitude');
title('Modulated BPSK Signal');
grid on;

# Spectrum domain
figure(4);
plot(f, (spectrum.^2));
xlabel('Frequency');
ylabel('Magnitude');
title(' Spectrum domain of the Modulated BPSK Signal');
grid on;

# Reciever
Carrier=cos(2 * pi * fc * t);
k = sqrt(2 / Tb) * mod_signal.* Carrier;
M=fftshift(fft (k)) / N;
H=abs(f)<bit_rate;
rec= real(ifft(ifftshift(H.*M))*N);


%decision code
for i=1:length(bits);
  if  rec(((i-1)*n+i*n)/(2*bit_rate))>=0;
    rec(((i - 1) * n)/bit_rate + 1:(i * n)/bit_rate) = 1;
   else
      rec(((i - 1) * n)/bit_rate + 1:(i * n)/bit_rate) = -1;
  end
end


%compare 1
error_bits_1=0;
 for i=1:length(bits);
  if  rec(((i-1)*n+i*n)/(2*bit_rate))!=sig(((i-1)*n+i*n)/(2*bit_rate));
    error_bits_1++;
  end
end
BER1=error_bits_1/length(bits);
disp(BER1);
