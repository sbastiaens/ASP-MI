%% ASP&MI Part 1

%% 1.1

% Generation of a signal
Fs = 200;
t = 0:1/Fs:1-1/Fs;
y = sin(2*pi*20*t)+ sin(2*pi*10*t);

% PSD defintion 1
Rxx = xcorr(y,'biased');
Rxxdft = abs(fftshift(fft(Rxx)));
freq = -Fs/2:Fs/length(Rxx):Fs/2-(Fs/length(Rxx));

figure(1)
plot(freq,Rxxdft);
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('PSD as the DTFT of the ACF (Definition 1)')

% PSD defintion 2
figure(2)
N = length(y);
xdft = fftshift(fft(y));
psdx = (1/(Fs*N)) * abs(xdft).^2;
freq = -Fs/2:Fs/length(y):Fs/2-(Fs/length(y));

plot(freq,(psdx))
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('PSD based on Definition 2')

%% 1.2

%a) Sunspot Time Series

load sunspot.dat
y_sun = sunspot(:,2);

[pxx,f] = periodogram(y_sun,[],[],1);    % Periodogram applied to sunspot
figure(3)
subplot(1,2,1)
plot(f,10*log10(pxx))
ylim([-20 60])
xlabel('Frequency (Cycles/Year)')
ylabel('Amplitude (dB)')
title('Periodogram of Sunspot Time Series')

%Same analysis for mean centered sunspot data:
y_mean = sunspot(:,2)-mean(sunspot(:,2));         % Remove mean
y_det = detrend(y_mean)+mean(sunspot(:,2));       % Added mean to see trend effect

[pxx_mean,f_mean] = periodogram(y_mean,[],[],1);
[pxx_trend,f_trend] = periodogram(y_det,[],[],1);

hold on 
plot(f_mean,10*log10(pxx_mean))
hold on
plot(f_trend, 10*log10(pxx_trend))
legend('Original', 'Without Mean', 'Without Trend ')

% Applying log to data
y_log = log10(y_sun);

k1 = find(y_log==-Inf);             % Remove variables equal to -Inf
for i=1:length(k1)
    y_log(k1(i)) = 0;
end
y_logmean = y_log-mean(y_log);
[pxx_log,f] = periodogram(y_log,[],[],1);

subplot(1,2,2)
plot(f,10*log10(pxx))
[pxx_logmean,f] = periodogram(y_logmean,[],[],1);
hold on
plot(f,10*log10(pxx_logmean))
ylim([-60 60])
xlabel('Frequency (Cycles/Year)')
ylabel('Amplitude (dB)')
title('Periodogram of Sunspot Time Series')
legend('Original','Log data without mean')

%b) EEG experiment

load('EEG_Data_Assignment1.mat')
clear pxx
clear f
[pxx, f] = periodogram(POz,[],5*1200,1200);    %Standard periogram

%For average periodogram window length L equal to (96000*L/80)samples

%Averaged Periodogram window length = 10s
var = 1;
for i=1:8
   x_10(i,:) = POz(var:var+11999);     
   var = var+12000;
end

for i=1:8
    [pxx_10(i,:), f_10(i,:)] = periodogram(x_10(i,:),[],5*1200,1200);
end
Average_10 = mean(pxx_10);
figure(4)
plot(f,10*log10(pxx)) 
hold on
plot(f,10*log10(Average_10))
xlim([0 100])
ylim([-160 -90])
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')
legend('Standard Periodogram','Averaged Periodogram window length 10s')

%Averaged Periodogram window length = 5s
var=1;
for i=1:16
   x_5(i,:) = POz(var:var+5999);
   var = var+6000;
end
for i=1:16
    [pxx_5(i,:), f_5(i,:)] = periodogram(x_5(i,:),[],5*1200,1200);
end
Average_5 = mean(pxx_5);

%Averaged Periodogram window length = 1s
var=1;
for i=1:80
   x_1(i,:) = POz(var:var+1199);
   var = var+1200;
end
for i=1:80
    [pxx_1(i,:), f_1(i,:)] = periodogram(x_1(i,:),[],5*1200,1200);
end
Average_1 = mean(pxx_1);

figure (5) 
plot(f,10*log10(Average_10))
hold on 
plot(f,10*log10(Average_5))
hold on 
plot(f,10*log10(Average_1))
xlim([0 100])
ylim([-140 -90])
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')
legend('Window length = 10s','Window length = 5s','Window length = 1s')

%% 1.3

% a) Validation of ACF and correlogram
y1= randn(1000,1);
t = (0:0.001:1);
fs = 1000;
y2 = 2*sin(2*pi*50*t)+ 0.5*sin(2*pi*100*t) + randn(size(t));
y2 = y2(1:end-1);
y3=filter(ones(9,1),[1],y1);

[r_bias, r_unbias, lag] = autocorr_scipt(y1);
[r_bias2, r_unbias2, lag2] = autocorr_scipt(y2);
[r_bias3, r_unbias3, lag3] = autocorr_scipt(y3);

pxx_bias = real(fftshift(fft(ifftshift(r_bias))));
pxx_unbias = real(fftshift(fft(ifftshift(r_unbias))));
pxx_bias2 = real(fftshift(fft(ifftshift(r_bias2))));
pxx_unbias2 = real(fftshift(fft(ifftshift(r_unbias2))));
pxx_bias3 = real(fftshift(fft(ifftshift(r_bias3))));
pxx_unbias3 = real(fftshift(fft(ifftshift(r_unbias3))));

f = linspace(-fs/2, fs/2, length(pxx_bias));
figure(6)
subplot(3,2,1)
plot(lag, r_unbias)
hold on
plot(lag, r_bias)
xlabel('Lags')
ylabel('Autocorrelation')
title('Autocorrelation of WGN')
legend('Unbiased Estimator', 'Biased Estimator')
subplot(3,2,3)
plot(lag2, r_unbias2)
hold on
plot(lag2, r_bias2)
xlabel('Lags')
ylabel('Autocorrelation')
title('Autocorrelation of noisy sinusoidal signal')
legend('Unbiased Estimator', 'Biased Estimator')
subplot(3,2,5)
plot(lag3, r_unbias3)
hold on
plot(lag3, r_bias3)
xlabel('Lags')
ylabel('Autocorrelation')
title('Autocorrelation of filtered WGN')
legend('Unbiased Estimator', 'Biased Estimator')

subplot(3,2,2)
plot(f, pxx_unbias)
hold on
plot(f, pxx_bias)
xlabel('Frequency (Hz)')
ylabel('PSD')
title('Correlogram spectral estimator of WGN')
legend('Unbiased Estimator', 'Biased Estimator')
subplot(3,2,4)
plot(f, pxx_unbias2)
hold on
plot(f, pxx_bias2)
xlabel('Frequency (Hz)')
ylabel('PSD')
title('Correlogram spectral estimator of noisy sinusoidal')
legend('Unbiased Estimator', 'Biased Estimator')
subplot(3,2,6)
plot(f, pxx_unbias3)
hold on
plot(f, pxx_bias3)
xlabel('Frequency (Hz)')
ylabel('PSD')
title('Correlogram spectral estimator of  filtered WGN')
legend('Unbiased Estimator', 'Biased Estimator')

%b) PSD estimate of different realisations

figure(7)
subplot(2,1,1)
t = (0:0.001:1);
t =t(1:end-1);
fs = linspace(0,500,1000);
for n=1:100
    w = randn(100,1)';
   x = 1*sin(2*pi*t*50) + 0.5*sin(2*pi*t*100)+5*randn(1000,1)';
    [r_bias, r_unbias, lag] = autocorr_scipt(x);
    pxx_bias = real(fftshift(fft(ifftshift(r_bias))));
    pxx_bias = pxx_bias(end/2:end);
    pxx_biasmatrix(n,:) = pxx_bias;
    hold on
    plot(fs, pxx_bias, 'b');
end
mn = mean(pxx_biasmatrix);
hold on
plot(fs,mn, 'r');
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('PSD estimates (different realisations and mean)')
subplot(2,1,2)
st = std(pxx_biasmatrix);
plot(fs,st);
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Standard deviation of the PSD estimate')

%c) PSD estimates in dB

figure (8)
subplot(2,1,1)
for n=1:100
    hold on
    plot(fs, 10*log10(pxx_biasmatrix(n,:)), 'b');
end
hold on
m = mean(10*log10(pxx_biasmatrix));
plot(fs, m, 'r')
xlabel('Frequency (Hz)')
ylabel('dB')
title('PSD estimates (different realisations and mean)')
subplot(2,1,2)
s = std(10*log10(pxx_biasmatrix));
plot(fs, s)
xlabel('Frequency (Hz)')
ylabel('dB')
title('Standard deviation of the PSD estimate')

%d) Verification of periodogram with more samples

n=0:50;
noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
x=exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+noise;
[pxx, f] = periodogram(x,rectwin(length(x)),128,1);
figure(9)
plot(f, 10*log10(pxx))
xlabel('Frequency (Hz)')
ylabel('Power/frequency (dB/Hz)')
title('Periodogram Power Spectral Density Estimate')

%e) MUSIC algorithm

figure(10)
[X,R] = corrmtx(x,14,'mod');
[S,F] = pmusic(R,2,[ ],1,'corr');
plot(F,S,'linewidth',2); set(gca,'xlim',[0.25 0.40]);
grid on; xlabel('Hz'); ylabel('Pseudospectrum');

%% 1.4

%b) PSD estimate of different model order compared to true PSD

x(1) = randn; x(2) = randn; x(3) = randn; x(4) = randn;
for n=5:10000
   x(n)=2.76*x(n-1)-3.81*x(n-2)+2.65*x(n-3)-0.92*x(n-4)+randn;
end
x=x(500:end);

%Ideal model
a = [1 -2.76 3.81 -2.65 0.92];
[H,F] = freqz(1,a, length(x));

%Model with estimated coefficients 
order = 2:14;
for i=1:length(order)
    [a_est, w(1, i)] = aryule(x, order(i));  %use Yule walker to determine coeff
    [res(:, i), F] = freqz(w(1, i)^(1/2), a_est, length(x)); %PSD estimate
end

figure(11) 
 v = [1 3 7 13];
 for i=1:4
     subplot(2,2,i)
     plot(F/pi,20*log10(abs(H)))
     hold on
     plot(F/pi,20*log10(abs(res(:,v(i)))))
     xlabel('Normalized Frequency')
     ylabel('PSD (dB/Hz)')
     legend('true PSD', 'estimated PSD')
     title(['PSD of ideal model and estimated model AR(',num2str(order(v(i))),')'])
 end

%% 1.5
% load('RAW.mat')
% Trial_1 = data(find(time==74.43):find(time==237.1));
% Trial_2 = data(find(time==245.2):find(time==489.5));
% Trial_3 = data(find(time==494.3):find(time==736.9));
% [RR1_1, RRIf_1] = ECG_to_RRI(Trial_1,1000);
% [RR1_2, RRIf_2] = ECG_to_RRI(Trial_2,1000);
% [RR1_3, RRIf_3] = ECG_to_RRI(Trial_3,1000);
% save('RRI_1.mat','RRI_1','RRIf_1');
% save('RRI_2.mat','RRI_2','RRIf_2');
% save('RRI_3.mat','RRI_3','RRIf_3');
load('RRI_1')
load('RRI_2')
load('RRI_3')
% t=1/RRIf_1:1/RRIf_1:length(RR1_1)/RRIf_1;
% figure
% plot(t,RR1_1)
% load('RRI_3')
% t=1/RRIf_2:1/RRIf_2:length(RR1_2)/RRIf_2;
% figure
% plot(t,RR1_2)
% t=1/RRIf_3:1/RRIf_3:length(RR1_3)/RRIf_3;
% figure
% plot(t,RR1_3)


%a) Standard and Averaged periodogram to three RRI trials

%Standard periodogram
[pxx1,f1] = periodogram(RR1_1,[],4*400,4);
[pxx2,f2] = periodogram(RR1_2,[],4*400,4);
[pxx3,f3] = periodogram(RR1_3,[],4*400,4);
figure(12)
subplot(3,2,1)
plot(f1,10*log10(pxx1))
xlabel('Frequency (Hz)')
ylabel('PSD')
title ('Standard Periodogram of Trial 1')
subplot(3,2,3)
plot(f2,10*log10(pxx2))
xlabel('Frequency (Hz)')
ylabel('PSD')
title ('Standard Periodogram of Trial 1')
subplot(3,2,5)
plot(f3,10*log10(pxx3))
xlabel('Frequency (Hz)')
ylabel('PSD')
title ('Standard Periodogram of Trial 1')

%Average periodogram
Window = [50 100 150];
subplot(3,2,2)
plot(f1,10*log10(pxx1))

for i=1:3
    [pxx_aver, w] = pwelch(RR1_1,hanning(Window(i)),0,4*400, 4);
    hold on
    plot(w,10*log10(pxx_aver))
end
xlabel('Frequency (Hz)')
ylabel('PSD')
title ('Trial 1 Averaged periodogram for different window length')
legend('Standard', 'Window length = 50', 'Window length = 100', 'Window length = 100')
subplot(3,2,4)
plot(f2,10*log10(pxx2))

for i=1:3
    [pxx_aver2, w2] = pwelch(RR1_2,hanning(Window(i)),0,4*400, 4);
    hold on
    plot(w2,10*log10(pxx_aver2))
end
xlabel('Frequency (Hz)')
ylabel('PSD')
title ('Trial 2 Averaged periodogram for different window length')
legend('Standard', 'Window length = 50', 'Window length = 100', 'Window length = 100')

subplot(3,2,6)
plot(f3,10*log10(pxx3))

for i=1:3
    [pxx_aver3, w3] = pwelch(RR1_3,hanning(Window(i)),0,4*400, 4);
    hold on
    plot(w3,10*log10(pxx_aver3))
end
xlabel('Frequency (Hz)')
ylabel('PSD')
title ('Trial 3 Averaged periodogram for different window length')
legend('Standard', 'Window length = 50', 'Window length = 100', 'Window length = 100')

%c)
RRtest = RR1_1.*hann(length(RR1_1))';
figure(13)
subplot(3,2,1)
plot(f1,10*log10(pxx1))
hold on
for i=1:30
    [d1_1, p1_1] = aryule(RR1_1, i);  %use Yule walker to determine coeff
    [H1_1,w_1] = freqz(sqrt(p1_1),d1_1, length(RR1_1), 4);
     plot(w_1,20*log10(abs(H1_1)))
     hold on
end
xlabel('Frequency (Hz)')
ylabel('AR spectrum')
title('AR spectrum of original data with standard periodogram ')
subplot(3,2,2)
plot(f1,10*log10(pxx1))
hold on
for i=1:20
    [d1_1, p1_1] = aryule(RRtest, i);  %use Yule walker to determine coeff
    [H1_1,w_1] = freqz(sqrt(p1_1),d1_1, length(RR1_1), 4);
    plot(w_1,20*log10(abs(H1_1)))
    hold on
end

xlabel('Frequency (Hz)')
ylabel('AR spectrum')
title('AR spectrum of windowed data with standard periodogram ')

RRtest2 = RR1_2.*hann(length(RR1_2))';
subplot(3,2,3)
plot(f2,10*log10(pxx2))
hold on
for i=1:32
    [d1_2, p1_2] = aryule(RR1_2, i);  %use Yule walker to determine coeff
    [H1_2,w_2] = freqz(sqrt(p1_2),d1_2, length(RR1_2), 4);
     plot(w_2,20*log10(abs(H1_2)))
     hold on
end
xlabel('Frequency (Hz)')
ylabel('AR spectrum')
title('AR spectrum of original data with standard periodogram ')
subplot(3,2,4)
plot(f2,10*log10(pxx2))
hold on
for i=1:20
    [d1_2, p1_2] = aryule(RRtest2, i);  %use Yule walker to determine coeff
    [H1_2,w_2] = freqz(sqrt(p1_2),d1_2, length(RR1_2), 4);
    plot(w_2,20*log10(abs(H1_2)))
    hold on
end
xlabel('Frequency (Hz)')
ylabel('AR spectrum')
title('AR spectrum of windowed data with standard periodogram ')

RRtest3 = RR1_3.*hann(length(RR1_3))';
subplot(3,2,5)
plot(f2,10*log10(pxx3))
hold on
for i=1:25
    [d1_3, p1_3] = aryule(RR1_3, i);  %use Yule walker to determine coeff
    [H1_3,w1_3] = freqz(sqrt(p1_3),d1_3, length(RR1_3), 4);
     plot(w1_3,20*log10(abs(H1_3)))
     hold on
end
xlabel('Frequency (Hz)')
ylabel('AR spectrum')
title('AR spectrum of original data with standard periodogram ')

subplot(3,2,6)
plot(f3,10*log10(pxx3))
hold on
for i=1:20
    [d1_3, p1_3] = aryule(RRtest3, i);  %use Yule walker to determine coeff
    [H1_3,w1_3] = freqz(sqrt(p1_3),d1_3, length(RR1_3), 4);
     plot(w1_3,20*log10(abs(H1_3)))
     hold on
end
xlabel('Frequency (Hz)')
ylabel('AR spectrum')
title('AR spectrum of windowed data with standard periodogram ')


    
%% 1.6 Robust Regression

%a) Singular Values and Square error
load('PCAPCR.mat')
s_x = svd(X);
s_noise = svd(Xnoise);
square_error = abs(s_x-s_noise).^2;

figure(14)
subplot(1,2,1)
stem(1:length(s_x), s_x)
hold on 
stem(1:length(s_x), s_noise)
xlabel('Singular Value Index')
ylabel('Singular Value Amplitude')
title('Singular Values of X and Xnoise')
legend('svd of X', 'svd of Xnoise')
subplot(1,2,2)
stem(1:length(square_error), square_error)
xlabel('Singular Value Index')
ylabel('Square Error')
title('Square error between each singular value of X and Xnoise')

%b) Low-rank approximation
[U S V] = svd(Xnoise);
k = 3;
X_noiseapprox = U(:,1:k)*S(1:k,1:k)*V(:,1:k)'; % rank-3 approximation
    
for i=1:10
   err(:,i) = abs(X(:,i) - Xnoise(:,i)).^2;
    err_mean(i) = mean(err(:,i));
     errapprox(:,i) = abs(X(:,i) - X_noiseapprox(:,i)).^2;
    err_meanapprox(i) = mean(errapprox(:,i));
end
figure(15)
stem(1:10,err_meanapprox)
hold on
stem(1:10,err_mean)
xlabel('Variables Index')
ylabel('Error')
title('Difference between the variables')
legend('X_{denoised}', 'X_{noise}')
diff1 = norm(X-X_noiseapprox,'fro');  
diff2 = norm(X-Xnoise,'fro');

%for different k
for k=1:10
    X_approx = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
    diff_k(k) = norm(X-X_approx,'fro'); 
    clear X_approx
end

figure(16)
plot(1:10, diff_k)
xlabel('Numer of components in the low-rank approximation') 
ylabel('Error between X and X_{denoised} (Frobenius norm)')

%c) OLS and PCR
k=3;
X_noiseapprox = U(:,1:k)*S(1:k,1:k)*V(:,1:k)'; 
B_ols = inv((Xnoise'*Xnoise)) * Xnoise'*Y;
B_pcr = V(:,1:k)*inv(S(1:k,1:k))*U(:,1:k)'*Y;

Y_ols = Xnoise*B_ols;
Y_pcr = X_noiseapprox*B_pcr;
Error_Yols = norm(Y-Y_ols,'fro');
Error_Ypcr = norm(Y-Y_pcr,'fro');

for t=1:10
    X_n = U(:,1:t)*S(1:t,1:t)*V(:,1:t)'; 
    B_p = V(:,1:t)*inv(S(1:t,1:t))*U(:,1:t)'*Y;
    Y_p = X_noiseapprox*B_p;
    Error_Yp(t) = norm(Y-Y_p,'fro');
    Error_Yo(t) = norm(Y-Y_ols,'fro');
    clear X_n
    clear B_p
    clear Y_p
end
figure(17)
subplot(1,2,1)
scatter(1:10, Error_Yp, 'x', 'LineWidth', 2);
hold on 
scatter(1:10, Error_Yo, 'x', 'LineWidth', 2);
xlabel('Number of components in low-rank Approximation')
ylabel('Estimation Error')
legend('PCR','OLS')
title('Estimation error between Y and the solutions found with PCR or OLS methods')
%With testing data
clear U
clear S
clear V
[U, S, V] = svd(Xtest);
k=3;
X_testapprox = U(:,1:k)*S(1:k,1:k)*V(:,1:k)'; 

Y_olstest = Xtest*B_ols;
Y_pcrtest = X_testapprox*B_pcr;
Error_Yolstest = norm(Ytest-Y_olstest,'fro');
Error_Ypcrtest = norm(Ytest-Y_pcrtest,'fro');

for t=1:10
    X_ntest = U(:,1:t)*S(1:t,1:t)*V(:,1:t)'; 
    B_ptest = V(:,1:t)*inv(S(1:t,1:t))*U(:,1:t)'*Ytest;
    Y_ptest = X_testapprox*B_ptest;
    Error_Yptest(t) = norm(Ytest-Y_ptest,'fro');
    Error_Yotest(t) = norm(Ytest-Y_olstest,'fro');
    clear X_ntest
    clear B_ptest
    clear Y_ptest
end
subplot(1,2,2)
scatter(1:10, Error_Yptest, 'x', 'LineWidth', 2);
hold on 
scatter(1:10, Error_Yotest, 'x', 'LineWidth', 2);
xlabel('Number of components in low-rank Approximation')
ylabel('Estimation Error')
legend('PCR','OLS')
title('Estimation error between Y_{test} and the solutions found with PCR or OLS methods')


%d) Efficacy assessment on test data
[OLSY_est, OLSYdata] = regval(B_ols);
[PCRY_est, PCRYdata] = regval(B_pcr);
err_ols = immse(OLSYdata,OLSY_est);
err_pcr = immse(PCRYdata,PCRY_est);