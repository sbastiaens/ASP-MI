%% ASP&MI Part 3 

clear all 
close all

%% 3.1

%a) ACLMS and CLMS

b1 = 1.5 + 1i;
b2 = 2.5 - 0.5i;
var = 1;
N = 1000;
for k = 1:100
    x(k,:) = sqrt(var/ 2) * (randn(1, N) + (1i * randn(1, N)));
end
y = zeros(100,N);
for k=1:100
    y(k,1) = x(k,1);
    for n = 2:N
        y(k,n) = x(k,n) + b1*x(k,n-1) + b2*conj(x(k,n-1));
    end
end
figure(1)
subplot(2,1,1)
scatter(real(x(1,:)), imag(x(1,:)), 30, "filled")
xlabel('Real part')
ylabel('Imaginary Part')
title('Plot of x(n)')
subplot(2,1,2)
scatter(real(y(1,:)), imag(y(1,:)), 30, "filled")
xlabel('Real part')
ylabel('Imaginary Part')
title('Plot of y(n)')

yCLMS =  zeros(100,N);
hCLMS = zeros(2,N+1,100);
eCLMS = zeros(100,N);
yACLMS =  zeros(100,N);
hACLMS = zeros(2,N+1,100);
eACLMS = zeros(100,N);
mu=0.1;
for k=1:100
    [yCLMS(k,:), ~, eCLMS(k,:)]= CLMS(x(k,:),y(k,:),2,N,mu, 'Identification');
    [yACLMS(k,:), ~, ~, eACLMS(k,:)]= ACLMS(x(k,:),y(k,:),2,N,mu, 'Identification');
end
Avg_CLMSe = mean(abs(eCLMS.^2));
Avg_ACLMSe = mean(abs(eACLMS.^2));
figure(2)
plot(1:N, 10*log10(Avg_CLMSe))
hold on
plot(1:N, 10*log10(Avg_ACLMSe))
legend('CLMS', 'ACLMS')
title('Learning curve of CLMS and ACLMS')
xlabel('Time index')
ylabel('Error in dB')

%b) Wind regimes

load('high-wind.mat')
N = length(v_east);
v_high = zeros(1,N);
for n=1:N
    v_high(n) = v_east(n) + 1i*v_north(n);
end

load('low-wind.mat')
v_low = zeros(1,N);
for n=1:N
    v_low(n) = v_east(n) + 1i*v_north(n);
end

load('medium-wind.mat')
v_medium = zeros(1,N);
for n=1:N
    v_medium(n) = v_east(n) + 1i*v_north(n);
end

rho_low = coeff_circularity(v_low);
rho_medium = coeff_circularity(v_medium);
rho_high = coeff_circularity(v_high);

figure(3)
subplot(3,1,1)
scatter(real(v_low), imag(v_low), 30, "filled")
xlabel('Real part')
ylabel('Imaginary part')
title('Scatter diagram of low wind regime')
grid on, grid minor
subplot(3,1,2)
scatter(real(v_medium), imag(v_medium), 30, "filled")
xlabel('Real part')
ylabel('Imaginary part')
title('Scatter diagram of medium wind regime')
grid on, grid minor
subplot(3,1,3)
scatter(real(v_high), imag(v_high), 30, "filled")
xlabel('Real part')
ylabel('Imaginary part')
title('Scatter diagram of high wind regime')
grid on, grid minor

wind = [v_low; v_medium; v_high];
mu = [0.1, 0.01, 0.001];
M = reshape(1:25,1,25);
for p = 1:3
    loweCLMS = zeros(max(M),N);
    loweACLMS = zeros(max(M),N);
    for k = 1:length(M)
        [~, ~, loweCLMS(k,:)]= CLMS(wind(p,:),wind(p,:),M(k),N,mu(p), 'Prediction');
        [~, ~, ~, loweACLMS(k,:)]= ACLMS(wind(p,:),wind(p,:),M(k),N,mu(p), 'Prediction');
    end
    for k=1:max(M)
       Avg_eclms(p,k) =mean(abs(loweCLMS(k,k:end).^2) );
       Avg_eaclms(p,k) =mean(abs(loweACLMS(k,k:end).^2) );
    end
end

figure(4)
for p=1:3
    subplot(3,1,p)
    plot(1:max(M), 10*log10(Avg_eclms(p,:)))
    hold on
    plot(1:max(M), 10*log10(Avg_eaclms(p,:)))
    legend('CLMS', 'ACLMS')
    title('Learning curve of CLMS and ACLMS')
    xlabel('Time index')
    ylabel('Error in dB')
    grid on, grid minor
end

%c) Complex Voltages: Balanced and Unbalanced System

fo = 50;
fs = 5000;
N = 1500;
n = 1:N;
V = ones(3,1);
phi = [0; -(2*pi/3); 2*pi/3];
delta = zeros(3,1);
v = V.*cos(2*pi*n*(fo/fs)+delta+phi);
clarkematrix = (sqrt(2/3))*[sqrt(2)/2, sqrt(2)/2, sqrt(2)/2; ...
                1, -1/2, -1/2;...
                0, sqrt(3)/2, -sqrt(3)/2];
vclarke = clarkematrix*v;
vfinal = vclarke(2,:) + vclarke(3,:)*1i;

figure(5)
subplot(1,3,1)
plot(real(vfinal), imag(vfinal))
xlabel('Real Part')
ylabel('Imaginary Part')
title('Balanced System with V = 1')
res1 = coeff_circularity(vfinal);

delta_unbal = [0; 0.2; 0.4];
Vc = [0.2, 0.4, 0.8, 3, 4];
subplot(1,3,2)
for q = 1:length(Vc)
    V_unbal = [1; 2; Vc(q)];
    v_unbal = V_unbal.*cos(2*pi*n*(fo/fs)+delta_unbal+phi);
    vclarke_unbal = clarkematrix*v_unbal;
    vfinal_unbal = vclarke_unbal(2,:) + vclarke_unbal(3,:)*1i;
    hold on
    plot(real(vfinal_unbal), imag(vfinal_unbal))
end
xlabel('Real Part')
ylabel('Imaginary Part')
legend(['V_c = ', num2str(Vc(1))], ['V_c  = ', num2str(Vc(2))], ['V_c  = ', num2str(Vc(3))], ['V_c  = ', num2str(Vc(4))])
title('Unbalanced System with varying V_c: V_b = 2, V_a = 1, \Delta_b = 0.2 and \Delta_c = 0.4')

delta_c = [0.1, 0.4, 0.8, 0.9];
subplot(1,3,3)
for q = 1:length(delta_c)
    delta_unbal = [0; 0.5; delta_c(q)];
    V_unbal = [1; 2; 4];
    v_unbal = V_unbal.*cos(2*pi*n*(fo/fs)+delta_unbal+phi);
    vclarke_unbal = clarkematrix*v_unbal;
    vfinal_unbal = vclarke_unbal(2,:) + vclarke_unbal(3,:)*1i;
    hold on
    plot(real(vfinal_unbal), imag(vfinal_unbal))
end
xlabel('Real Part')
ylabel('Imaginary Part')
legend(['\Delta_c = ', num2str(delta_c(1))], ['\Delta_c = ', num2str(delta_c(2))], ['\Delta_c = ', num2str(delta_c(3))], ['\Delta_c = ', num2str(delta_c(4))])
title('Unbalanced System with varying \Delta_c: \Delta_b = 0.5 and V_a = 1, V_b = 2 and V_c = 4')

res = coeff_circularity(vfinal_unbal);

%e) Frequency estimation with CLMS and ACLMS

mu = 0.05;
M = 1;
[~, h1, e]= CLMS(vfinal,vfinal,M,N,mu,'Prediction');
h1 = conj(h1);
fo_CLMS = (fs/(2*pi))*atan(imag(h1(2:end))./real(h1(2:end)));
figure(6)
subplot(1,2,1)
plot(1:N, 10*log10(abs(e.^2)))
[~, h, g, err]= ACLMS(vfinal,vfinal,M,N,mu,'Prediction');
h = conj(h);
g = conj(g);
fo_ACLMS = (fs/(2*pi))*atan(sqrt((imag(h(2:end)).^2 - (abs(g(2:end))).^2))./real(h(2:end)));
hold on 
plot(1:N, 10*log10(abs(err.^2)))
ylabel('Squared error prediction')
xlabel('Sample number')
legend('CLMS', 'ACLMS')
title('Prediction error of frequency estimation of the \alpha - \beta voltages')

subplot(1,2,2)
plot(1:N, fo_CLMS)
hold on
plot(1:N, fo_ACLMS)
xlabel('Sample number')
ylabel('Frequency Estimation')
title('Nominal frequency estimation')
legend('CLMS','ACLMS')


[~, h1, e]= CLMS(vfinal_unbal,vfinal_unbal,M,N,mu,'Prediction');
h1 = conj(h1);
fo_CLMSunbal = (fs/(2*pi))*atan(imag(h1(2:end))./real(h1(2:end)));
figure(7)
plot(1:N, 10*log10(abs(e.^2)))
[~, h, g, err]= ACLMS(vfinal_unbal,vfinal_unbal,M,N,mu,'Prediction');
h = conj(h);
g = conj(g);
fo_ACLMSunbal = (fs/(2*pi))*atan(sqrt((imag(h(2:end)).^2 - (abs(g(2:end))).^2))./real(h(2:end)));
hold on
plot(1:N, 10*log10(abs(err.^2)))
ylabel('Squared error prediction')
xlabel('Sample number')
legend('CLMS', 'ACLMS')
title('Prediction error of frequency estimation of the \alpha - \beta voltages')

figure(8)
plot(1:N, fo_CLMSunbal)
hold on
plot(1:N, fo_ACLMSunbal)
ylim([-160 160])
xlabel('Sample number')
ylabel('Frequency')
title('Nominal frequency estimation')
legend('CLMS','ACLMS')

delta_unbal = [0; 0.2; 0.4];
V_unbal = [1; 0.5; 1.5];
v_unbal = V_unbal.*cos(2*pi*n*(fo/fs)+delta_unbal+phi);
vclarke_unbal = clarkematrix*v_unbal;
vfinal_unbal(q,:) = vclarke_unbal(2,:) + vclarke_unbal(3,:)*1i;
[~, h1, e]= CLMS(vfinal_unbal(q,:),vfinal_unbal(q,:),M,N,mu,'Prediction');
h1 = conj(h1);
fo_CLMSunbal = (fs/(2*pi))*atan(imag(h1(2:end))./real(h1(2:end)));

figure(9)
subplot(1,2,1)
plot(1:N, 10*log10(abs(e.^2)))
[~, h, g, err]= ACLMS(vfinal_unbal(q,:),vfinal_unbal(q,:),M,N,mu,'Prediction');
h = conj(h);
g = conj(g);
fo_ACLMSunbal = (fs/(2*pi))*atan(sqrt((imag(h(2:end)).^2 - (abs(g(2:end))).^2))./real(h(2:end)));
hold on
plot(1:N, 10*log10(abs(err.^2)))
ylabel('Squared error prediction')
xlabel('Sample number')
legend('CLMS', 'ACLMS')
title('Prediction error of frequency estimation of the \alpha - \beta voltages')
subplot(1,2,2)
plot(1:N, fo_CLMSunbal)
hold on
plot(1:N, fo_ACLMSunbal)
ylim([-160 160])
xlabel('Sample number')
ylabel('Frequency')
title('Nominal frequency estimation')
legend('CLMS','ACLMS')
ylim([0, 120])

    
%% 3.2

clear h
clear w
clear y 

% a) FM signal with AR power spectrum
 
var = 0.05; 
N = 1500;
noise = sqrt(var/ 2) * (randn(1, N) + (1i * randn(1, N)));

fs = 1000;
for n = 1:N
   if n <=500
    f(n) = 100;
    phi(n) =  100*n;
   elseif n >= 1001
     f(n) = 100 + ((n-1000)/25)^(2);
     phi(n) = ((1/3)*(n^(3)/(25^(2)))) - ((1.6)*n^(2)) + (1700*n);  
   else
       f(n) = 100 + ((n-500)/2);
     phi(n) = (-150*n) + (0.25*n^(2));
   end
end
y = exp(1i*((2*pi)/fs)*phi')+noise;

figure(10)
order = [1, 10];
 for p = 1:2
     a = aryule(y, order(p)); 
     [h, w] = freqz(1, a, N, fs);
     subplot(1,2,p)
     plot(w,20*log10(abs(h)))
     xlabel('Frequency (Hz)')
     ylabel('Power Spectral Density (PSD) dB')
      title(['AR(',num2str(order(p)),') model of y'])
      grid on, grid minor
 end
 
 figure(11)
  order = [1, 10];
  t = 1;
 for p = 1:2
     a = aryule(y(1:500), order(p)); 
     [h, w] = freqz(1, a, 500, fs);
     subplot(2,3,t)
     plot(w,20*log10(abs(h)))
     xlabel('Frequency (Hz)')
     ylabel('Power Spectral Density (PSD) dB')
     title(['AR(',num2str(order(p)),') model of y(1:500)'])
     grid on, grid minor
     a = aryule(y(501:1000), order(p)); 
     [h, w] = freqz(1, a, 500, fs);
     subplot(2,3,t+1)
     plot(w,20*log10(abs(h)))
     xlabel('Frequency (Hz)')
     ylabel('Power Spectral Density (PSD) dB')
     title(['AR(',num2str(order(p)),') model of y(501:1000)'])
     grid on, grid minor
     a = aryule(y(1001:1500), order(p)); 
     [h, w] = freqz(1, a, 500, fs);
     subplot(2,3,t+2)
     plot(w,20*log10(abs(h)))
     xlabel('Frequency (Hz)')
     ylabel('Power Spectral Density (PSD) dB')
     title(['AR(',num2str(order(p)),') model of y(1001:1500)'])
     grid on, grid minor
     t = 4;
 end
 
 %b) CLMS to estimate AR coefficients

clear h
clear w

figure(12)
subplot(2,3,[1,3])
plot(1:N,f)
xlabel('Time index')
ylabel('Frequency (Hz)')
grid on, grid minor
title('Original time-frequency spectrum')

M = 1;
mu = [0.1, 0.05, 0.01];
for b = 1:length(mu) 
    for n = 1:N
     [~, h1, ~]= CLMS(y,y,M,N,mu(b),'Prediction');
    % Run complex-valued LMS algorithm to estimate AR coefficient a?1(n)
     [h ,w]= freqz(1 , [1; -conj(h1(n))], 1024, fs);  % Compute power spectrum
     Ht(:, n) = abs(h).^2; % Store it in a matrix;
    end
    % Remove outliers in the matrix H
    medianH = 50*median(median(Ht));
    Ht(Ht > medianH) = medianH;
    % Plot time-frequency diagram
    subplot(2,3,b+3)
    surf(1:N, w, Ht, "LineStyle", "none");
    view(2);
    xlabel('Time index (n)')
    ylabel('Frequency (Hz)')
    title(['Estimated time-frequency spectrum: \mu =',num2str(mu(b))])
end

%% 3.3

%c) DFT-CLMS algorithm

% FM signal from 3.2
clear x

var = 0.05; 
N = 1500;
    noise = sqrt(var/ 2) * (randn(1, N) + (1i * randn(1, N)));
fs = 1500;
for n = 1:N
   if n <=500
    f(n) = 100;
    phi(n) =  100*n;
   elseif n >= 1001
     f(n) = 100 + ((n-1000)/25)^(2);
     phi(n) = ((1/3)*(n^(3)/(25^(2)))) - ((1.6)*n^(2)) + (1700*n);  
   else
       f(n) = 100 + ((n-500)/2);
     phi(n) = (-150*n) + (0.25*n^(2));
   end
end

y = exp(1i*((2*pi)/fs)*phi')+noise;
M = 1;
mu = 0.05;
K = 2080;
F = zeros(N,K);

for n = 1:N
    for  p = 0:K-1
        F(n,p+1) = exp((1i*2*pi*n*p)/K);
    end
    x(:,n) = (1/K)*F(n,:)';
end
gamma = [0.5, 1, 3];

figure(13)
for i=1:length(gamma)
    [~, hDFT, ~]= DTFCLMS(x,y',K,N,mu,gamma(i));
    H = abs(hDFT).^2; % Store it in a matrix;
    w = (0:(K-1)) .* (fs / K);
    % Remove outliers in the matrix H
    medianH = 50*median(median(H));
    H(H > medianH) = medianH;
    % Plot time-frequency diagram
    subplot(1,3,i)
    surf(1:N, w, H, "LineStyle", "none");
    view(2);
    xlabel('Time index (n)')
    ylabel('Frequency (Hz)')
    ylim([0, 500])
    colorbar
    title(['Estimated time-frequency spectrum: \gamma =',num2str(gamma(i))])
end

 %d) DFT-CLMS for EEG signal
 
 load ('EEG_Data_Assignment1')
 a = 500;
 POz = POz(a:a+1200-1);
 N = length(POz);
 K = 2080;
 fs = 1200;
 F = zeros(N,K);
 mu = 0.1;
 for n = 1:N
    for  p = 0:K-1
        F(n,p+1) = exp((1i*2*pi*n*p)/K);
    end
    x(:,n) = (1/K)*F(n,:)';
 end
 
figure(14)

[~, hDFT, ~]= DTFCLMS(x,POz',K,N,mu,0);

H = abs(hDFT).^2; % Store it in a matrix;
w = (0:(K-1)) .* (fs / K);

% Remove outliers in the matrix H
medianH = 50*median(median(H));
H(H > medianH) = medianH;
% Plot time-frequency diagram
surf(1:N, w, H, "LineStyle", "none");
view(2);
xlabel('Time index (n)')
ylabel('Frequency (Hz)')
ylim([0, 80])
title('Time-frequency diagram with \mu = 0.1')
colorbar