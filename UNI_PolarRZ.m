total_BER2=[];
bits = randi([0, 1], 1, 10000);
n=800;
T=length(bits);
N=n*length(bits);
dt=T/N;
t=0:dt:(T-dt);
sig=zeros(1,length(t));
for i=1:length(bits);
  if bits(i)==1
    sig((i-1)*n+2:i*n-n/2)=1.2;

  end
end

#Plot Signal
figure(7);
plot(t,sig),axis([0,20,-2,2]);
xlabel("Time");
ylabel("Signal");
title("UNI-polar RZ");
grid on;

#EyeDiagram
eyediagram(sig,800,1,400);
xlabel("Time");
ylabel("Amplitude");
title("Unipolar RZ Eyediagram");

#spectral domain
df = 1/T;
fs= 1/dt;
if (rem(N,2)==0)
f=- (0.5 * fs) : df : (0.5*fs - df) ;
else
f=- (0.5*fs-0.5*df) : df : (0.5*fs -0.5*df) ;
end
M=fftshift(fft(sig))/N;
figure(9);

PSD_y=power(abs(M),2);
plot (f,PSD_y*1000),axis([-6,6,0,1]);
xlabel("Frequency");
ylabel("Power");
title("Unipolar RZ Spectral");

H=abs(f)<2;
k = real(ifft(ifftshift(H.*M))*N);

#Reciever without noise
#decision code 1
for i=1:length(bits);
  if  k(((i-1)*n+i*n-n/2)/2)>=0.6;
    k((i-1)*n+1:i*n-n/2)=1.2;
     k((i*n-n/2)+1:i*n)=0;
   else
      k((i-1)*n+1:i*n)=0;
  end
end


#compare 1
error_bits_1=0;
 for i=1:length(bits);
  if  k(((i-1)*n+i*n-n/2)/2)!=sig(((i-1)*n+i*n-n/2)/2);
    error_bits_1++;
  end
 end
BER1=error_bits_1/length(bits);
disp(BER1);

#Reciever with noise
for sigma = 0:0.12:1.2;
noise=sigma*randn(1,length(t));
k1=k+noise;

#decision code 2
 for i=1:length(bits);
  if  k1(((i-1)*n+i*n-n/2)/2)>=0.6;
    k1((i-1)*n+1:i*n-n/2)=1.2;
     k1((i*n-n/2)+1:i*n)=0;
   else
      k1((i-1)*n+1:i*n)=0;
  end
 end

#compare 2
error_bits_2=0;
 for i=1:length(bits);
  if  k1(((i-1)*n+i*n)/2)!=sig(((i-1)*n+i*n)/2);
   error_bits_2++;
  end
 end
BER2=error_bits_2/length(bits);
total_BER2 =[total_BER2 BER2];
disp(BER2);
end

figure(17)
sigma =[0:0.12:1.2];
semilogy(sigma,total_BER2);
hold on;

