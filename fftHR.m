function   [HR, freq, P1, f] = fftHR(x, Fs)

T = 1/Fs;            
N = length(x);         
t = (0:N-1)*T;      

seq=0:N-1;  

df=Fs/N; 
f=(0:1:N/2)*df;  

Y = fft(x); 
P2 = abs(Y/N); 
P1 = P2(1:N/2+1); 
P1(2:end-1) = P1(2:end-1)*2;  

[~, loca]=max(P1); 
freq = f(loca);
HR = roundn(freq*60,-2);

end

