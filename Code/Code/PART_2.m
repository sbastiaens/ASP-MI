%% ASP&MI Part 2 

%% 2.1
clear all 
close all

% b) LMS algorithm
mu = [0.05, 0.01];
noise = sqrt(0.25)*randn(1000,1);
x(1) = randn; x(2) = randn;
for n=3:1000
   x(n)=0.1*x(n-1)+0.8*x(n-2)+noise(n);
end
order = 2;
N=length(x); 

 figure(1)
 for t = 1:length(mu) 
    w = zeros(order,1);
    w_ad = zeros(order,N);
    LMS = zeros(1,N);
    n=1;
    clear xsum
    clear e
    e(1) = x(1);
     for i=order+1:N
       xsum = [x(i-1) x(i-2)];
       LMS(i) = w'*xsum'; 
       e(i) = x(i) - LMS(i);
       w = w + mu(t)*e(i)*xsum'; 
       w_ad(:,n) = w;
       n = n+1;
     end 
     hold on
     plot(1:1000, 10*log10(e.^2))
 end   
 xlabel('Sample number')
 ylabel('Squared error prediction (dB)')
 legend('\mu = 0.05', '\mu = 0.01')
 
%100 realisations
 noise = sqrt(0.25)*randn(100,1000);
for k=1:100
    for n=3:1000
       x(k,1) = randn; x(k,2) = randn;
       x(k,n)=0.1*x(k,n-1)+0.8*x(k,n-2)+noise(k,n);
    end
end

figure(2)
 for t = 1:length(mu) 
    clear e
    for k=1:100
    w = zeros(order,1);
    w_ad = zeros(order,N);
    LMS = zeros(1,N);
    n=1;
    clear xsum
         for i=order+1:N
           xsum = [x(k,i-1) x(k,i-2)];
           LMS(i) = w'*xsum'; 
           e(k,i) = x(k,i) - LMS(i);
           w = w + mu(t)*e(k,i)*xsum'; 
           w_ad(:,n) = w;
           n = n+1;
         end
    end
    Avg_e(t,:) = mean((e.^2));
    hold on
    plot(1:1000, (10*log(Avg_e(t,:))))
 end 
 xlabel('Sample number')
 ylabel('Squared error prediction (dB)')
 legend('\mu = 0.05', '\mu = 0.01')
 
%c) Misadjustement Calculations

%Estimated misadjustement
for t=1:length(mu)
    sig = 0.25;
    MSE_m(t)= mean(Avg_e(t,end-600:end));
    mis(t) = (MSE_m(t)-sig)/sig;  
end

%Theoretical LMS misadjustement
Rxx = (25/54)*[2,1; 1,2];
for t=1:length(mu)
    mis_LMS(t) = (mu(t)/2)*trace(Rxx);
end

%d) Steady state values of adaptive filter coeffficients

 for t = 1:length(mu) 
    clear e
    for k=1:100
    w = zeros(order,1);
    w_ad = zeros(order,N);
    LMS = zeros(1,N);
    n=1;
    clear xsum
         for i=order+1:N
           xsum = [x(k,i-1) x(k,i-2)];
           LMS(i) = w'*xsum'; 
           e(k,i) = x(k,i) - LMS(i);
           w = w + mu(t)*e(k,i)*xsum'; 
           w_ad(:,n) = w;
           n = n+1;
         end
       coeff1(k,:) = w_ad(1,:);
       coeff2(k,:) = w_ad(2,:);
    end
    Avg_coeff1(t,:) = mean(coeff1);
    Avg_coeff2(t,:) = mean(coeff2);
 end 
 
figure(3)
subplot(1,2,1)
plot(1:998, Avg_coeff1(1,1:end-2))
hold on
line('XData', [0 1000], 'YData', [0.8 0.8], 'LineStyle', '-. ','LineWidth', 2, 'Color','r')
plot(1:998, Avg_coeff2(1,1:end-2))
xlabel('Sample Number')
ylabel('Coefficient Value')
title('Values of weights with \mu = 0.05')
grid on, grid minor
line('XData', [0 1000], 'YData', [0.1 0.1], 'LineStyle', '-. ','LineWidth', 2, 'Color','b')
legend('w(1), estimate of a1','True a1','w(2), estimate of a2','True a2', 'Orientation', 'horizontal')
axis([0 1000 0 1]) 
subplot(1,2,2)
plot(1:998, Avg_coeff1(2,1:end-2))
hold on
line('XData', [0 1000], 'YData', [0.8 0.8], 'LineStyle', '-. ','LineWidth', 2, 'Color','r')
plot(1:998, Avg_coeff2(2,1:end-2))
xlabel('Sample Number')
ylabel('Coefficient Value')
title('Values of weights with \mu = 0.01')
grid on, grid minor
line('XData', [0 1000], 'YData', [0.1 0.1], 'LineStyle', '-. ','LineWidth', 2, 'Color','b') 
legend('w(1), estimate of a1','True a1','w(2), estimate of a2','True a2','Orientation', 'horizontal')
axis([0 1000 0 1]) 

%Percentage difference 
Perc_coeff1 = ((Avg_coeff1(:,500:end-2) - 0.1)/0.1) * 100; 
PercAvg1 = mean(Perc_coeff1,2);
Perc_coeff2 = ((Avg_coeff2(:,500:end-2) - 0.8)/0.8) * 100; 
PercAvg2 = mean(Perc_coeff2,2);

%e) Leaky LMS
gamma = [0.1, 0.5, 0.9];
count =1;
for t = 1:length(mu)
     for p = 1:length(gamma) 
        clear e
        for k=1:100
        w = zeros(order,1);
        w_ad = zeros(order,N);
        LMS = zeros(1,N);
        n=1;
        clear xsum
             for i=order+1:N
               xsum = [x(k,i-1) x(k,i-2)];
               LMS(i) = w'*xsum'; 
               e(k,i) = x(k,i) - LMS(i);
               w = (1 - mu(t)*gamma(p))*w + mu(t)*e(k,i)*xsum'; 
               w_ad(:,n) = w;
               n = n+1;
             end
           coeff1(k,:) = w_ad(1,:);
           coeff2(k,:) = w_ad(2,:);
        end
        Avg_coeff1(count,:) = mean(coeff1);
        Avg_coeff2(count,:) = mean(coeff2);
        count = count+1;
     end
end
 
 figure(4)
 for i = 1:6
 subplot(2,3,i)
 plot(1:998, Avg_coeff1(i,1:end-2))
 line('XData', [0 1000], 'YData', [0.8 0.8], 'LineStyle', '-. ','LineWidth', 2, 'Color','r')
hold on
 plot(1:998, Avg_coeff2(i,1:end-2))
line('XData', [0 1000], 'YData', [0.1 0.1], 'LineStyle', '-. ','LineWidth', 2, 'Color','b') 
legend('w(1), estimate of a1','True a1','w(2), estimate of a2','True a2')
xlabel('Sample Number')
ylabel('Coefficient Value')
grid on, grid minor
    if i<= 3
        title(['Values of weights with \mu = 0.05 and \gamma =', num2str(gamma(i))])
    else
        title(['Values of weights with \mu = 0.01 and \gamma =', num2str(gamma(i-3))])
    end
end

 
%% 2.2

clear x
clear noise

noise = sqrt(0.5)*randn(100,1000);
for k =1:100
    for n=2:1000
        x(k,1) = noise(k,1);
        x(k,n) = 0.9*noise(k,n-1)+noise(k,n);
    end
end

% a) GASS algorithm

N = length(x);
order = 1;
ro = 0.005;
method = {'Benveniste','Ang','Matthews'};
for t=1:5
    clear e
    for k=1:100
        clear mu
        w = zeros(order+1,1);
        LMS = zeros(1,N);
        if t<=3
          mu(1) = 0.1;
          mu(2)=0.1;
        elseif t==4
            mu(2) = 0.01;
        else
            mu(2) = 0.05;
        end
        term = 0;
        clear xsum
        n=1;
        e(k,1) = 0;
        xsum = zeros(k,2);
         for i=order+1:N
           xsum(k,:) = [noise(k,i), noise(k,i-1)];
           LMS(i) = w'*xsum(k,:)'; 
           e(k,i) = x(k,i) - LMS(i);
           w = w + mu(i)*e(k,i)*xsum(k,:)'; 
           w_ad(:,n) = w;
           n = n+1;
           if t<=3
               term = gassTerm(i,k,xsum(k,:),mu,e,term,method{t});
               mu(i+1) = mu(i) + ro*e(k,i)*xsum(k,1)'*term;
               save(t,i+1) = mu(i+1);
           else
               mu(i+1) = mu(i);
           end
         end
         coeff1(k,:) = w_ad(2,:);
    end 
     Avg_coeff1(t,:) = mean(coeff1);
    Avg_e(t,:) = mean((e.^2)); 
end
figure(5)
for i=1:5
    hold on
plot(1:1000,10*log10(Avg_e(i,:)))
end
ylabel('Squared Prediction Error')
xlabel('Sample Number')
title('Squared Prediction Error for GASS (\mu = 0.1) and LMS algoritm with fixed \mu')
figure(6)
 plot(1:1000, 0.9-(Avg_coeff1(1,:)))
 hold on
 plot(1:1000, 0.9-(Avg_coeff1(2,:)))
 hold on
 plot(1:1000, 0.9-(Avg_coeff1(3,:)))
 hold on
 plot(1:1000, 0.9-(Avg_coeff1(4,:)))
 hold on
 plot(1:1000, 0.9-(Avg_coeff1(5,:)))
 xlabel('Sample Number')
 ylabel('Weight error (0.9-w(n))')
 legend('Benveniste','Ang','Matthews','\mu = 0.01', '\mu = 0.05')
 title('Weight error curves for GASS (\mu = 0.1) and LMS algorithms')
 
 
%c) GNGD algorithm

clear x
clear noise
close all

noise = sqrt(0.5)*randn(100,1000);
for k =1:100
    for n=2:1000
        x(k,1) = noise(k,1);
        x(k,n) = 0.9*noise(k,n-1)+noise(k,n);
    end
end
N = length(x);
order = 1;
ro = 0.005;
beta=1;
for t=1:5
    clear e
    for k=1:100
        clear mu
        clear eta
        w = zeros(order+1,1);
        LMS = zeros(1,N);
          mu(1) = 0.05;
          mu(2)=0.05;
          eta(1) =10;
          eta(2) = 10;
        term = 0;
        clear xsum
        n=1;
        e(k,1) = 0;
        xsum = zeros(k,2);
         for i=order+1:N
           xsum(k,:) = [noise(k,i), noise(k,i-1)];
           LMS(i) = w'*xsum(k,:)'; 
           e(k,i) = x(k,i) - LMS(i);
           if t==1
            w = w + (beta/(eta(i)+(xsum(k,:)*xsum(k,:)')))*e(k,i)*xsum(k,:)'; 
            w_ad(:,n) = w;
            n = n+1;
            eta(i+1) = eta(i) - ro*mu(1)*((e(i)*e(i-1)*xsum(k,1)'*xsum(k,2))/((eta(i-1)+(xsum(k,1)'*xsum(k,2)))^2))*xsum(k,1);
            elseif t==2
                w = w + mu(i)*e(k,i)*xsum(k,:)'; 
                w_ad(:,n) = w;
                n = n+1;
                term = gassTerm(i,k,xsum,mu,e,term,'Benveniste');
                mu(i+1) = mu(i) + ro*e(k,i)*xsum(k,1)*term;

            end

         end
         coeff1(k,:) = w_ad(2,:);
    end 
     Avg_coeff1(t,:) = mean(coeff1);
    Avg_e(t,:) = mean((e.^2)); 
end

figure(7)
for i=1:2
    hold on
plot(1:1000,10*log10(Avg_e(i,:)))
end
ylabel('Squared Prediction Error')
xlabel('Sample Number')
title('Squared Prediction Error for Benveniste \mu = 0.05 and GNGD with \eta = 10')

figure(8)
 plot(1:1000, 0.9-(Avg_coeff1(1,:)))
 hold on
 plot(1:1000, 0.9-(Avg_coeff1(2,:)))
 xlabel('Sample Number')
 ylabel('Weight error (0.9-w(n))')
 legend('GNGD','Benveniste')
 
%% 2.3

clear x
clear noise
noise = zeros(100,1000);
x = zeros(100,1000);
wfr = 0.01*pi;
v = randn(100,1000);
for k=1:100
    for n=1:1000
       x(k,n) = sin(n*wfr);
    end
end
for k=1:100
    for n=3:1000
       noise(k,n) = v(k,n)+0.5*v(k,n-2);
    end
end
s = x + noise;
s_mean = mean(s);
M = [5, 10, 15, 20];
LMS_mean = zeros(125,1000);
delay = reshape (1:25,1,25);

%a) ALE

figure(8)
for i=1:4
    LMS = ALE(s,delay(i),M(1));
    subplot(2,2,i)
    p1 = plot(1:1000, s, "Color", "blue");
    hold on
    p2 = plot(1:1000, LMS, "Color", "red");
    hold on 
    p3 = plot(1:1000, x(1,:), "Color", "yellow");
    leg = legend([p1(1),p2(1),p3(1)],'s(n)','$\hat{x}(n)$','x(n)');
    set(leg,'Interpreter','latex');
    title(['ALE of 100 realisations with \Delta = ', num2str(delay(i)),' and M = 5'])
    ylabel('Amplitude')
    xlabel('Sample Number')
end

%b) Dependance between delay and MSPE + effect of M 

count=1;
for i = 1:length(M)
    for j = 1:length(delay)
        LMS = ALE(s,delay(j),M(i));
        e_mean(count,:) = mean((x-LMS).^2);
        count = count+1;
    end
end
figure(9)
subplot(1,2,1)
LMS = ALE(s,25,5);
p1 = plot(1:1000, s, "Color", "blue");
hold on
p2 = plot(1:1000, LMS, "Color", "red");
hold on 
p3 = plot(1:1000, x(1,:), "Color", "yellow");
leg = legend([p1(1),p2(1),p3(1)],'s(n)','$\hat{x}(n)$','x(n)');
set(leg,'Interpreter','latex');
title(['ALE of 100 realisations with \Delta = 25 and M = 5'])
ylabel('Amplitude')
xlabel('Sample Number')
subplot(1,2,2)
LMS = ALE(s,3,20);
p1 = plot(1:1000, s, "Color", "blue");
hold on
p2 = plot(1:1000, LMS, "Color", "red");
hold on 
p3 = plot(1:1000, x(1,:), "Color", "yellow");
leg = legend([p1(1),p2(1),p3(1)],'s(n)','$\hat{x}(n)$','x(n)');
set(leg,'Interpreter','latex');
title(['ALE of 100 realisations with \Delta = 3 and M = 20'])
ylabel('Amplitude')
xlabel('Sample Number')
    
for i=1:100
    MSPE(i)  = mean(e_mean(i,:));
end
figure(10)
count = 25;
for i=1:4
    plot(1:25,MSPE(1, count-24:count))
    hold on
    count = count +25;
end
legend('M = 5','M = 10','M = 15','M = 20')
title('MSPE as a function of \Delta for different M values')
ylabel('MSPE')
xlabel('Delay (\Delta)')

delay = 3;
M = reshape(1:20,1,20);
for i=1:length(M)
    LMS1 = ALE(s,delay,M(i));
    e_mean1(i,:) = mean((x-LMS1).^2);
    MSPE1(i)  = mean(e_mean1(i,:));
end

figure(11)
plot(1:20, MSPE1)
title('MSPE as a function of M for \Delta = 3 and \mu = 0.01')
ylabel('MSPE')
xlabel('Filter Length (M)')

M = [5, 10, 15, 20];

%c) ANC 

epsilon = noise * 0.9 + 0.05;

LMS = ANC(s,epsilon,5);
LMS_ALE = ALE(s, 3, 5);

figure(12)
subplot(1,2,2)
p1 = plot(1:1000, s, "Color", "blue");
hold on
p2 = plot(1:1000, LMS, "Color", "red");
hold on 
p3 = plot(1:1000, x(1,:), "Color", "yellow");
leg = legend([p1(1),p2(1),p3(1)],'s(n)','$\hat{x}(n)$','x(n)');
set(leg,'Interpreter','latex');
e_ANC = mean((x-LMS).^2);
MSPE_ANC  = mean(e_ANC);
title('ANC for 100 realisations with M = 5 and \mu = 0.01')
xlabel('Sample Number')
ylabel('Amplitude')

subplot(1,2,1)
p1 = plot(1:1000, s, "Color", "blue");
hold on
p2 = plot(1:1000, LMS_ALE, "Color", "red");
hold on 
p3 = plot(1:1000, x(1,:), "Color", "yellow");
leg = legend([p1(1),p2(1),p3(1)],'s(n)','$\hat{x}(n)$','x(n)');
set(leg,'Interpreter','latex');
e_ALE = mean((x-LMS_ALE).^2);
MSPE_ALE  = mean(e_ALE);
title('ALE for 100 realisations with \Delta = 3, M = 5 and \mu = 0.01')
xlabel('Sample Number')
ylabel('Amplitude')
count = 0;

for i=1:4
    MSPE_delay3(i) = MSPE(3+count);
    count = count +25;
end

for i=1:length(M)
    LMS = ANC(s,epsilon,M(i));
    e_ANC = mean((x-LMS).^2);
    MSPE_ANC(i)  = mean(e_ANC);
end


%d) EEG data Assignment 2
N=length(POz);
f=50;
ts=1/fs;
t = (1:N)*ts;
sinusoid=sin(2*pi*f*t);
 noise = randn(1,length(POz));
 ref_input =  sinusoid+noise;
M = [5, 18, 50];
mu = [0.01, 0.005, 0.001];

LMS = zeros(9,length(POz));
count =1;
for p = 1:length(M) 
    for d = 1:length(mu) 
        noise_est = zeros(1,length(POz));
        clear u
         w = zeros(M(p),1); 
        for n=M(p):length(POz)
            for i=1:M(p)
                u(i,n) = ref_input(n - i + 1);
            end
            noise_est(n) = w'*u(:,n);
            LMS(count,n) = POz(n) - noise_est(n);
            w = w + mu(d)*LMS(count,n)*u(:,n);
        end  
     count = count+1;
    end
end
L=5000;
noverlap =round(0.5*L);
K=16000;

figure(13)
subplot(1,2,1)
spectrogram(POz, rectwin(L), noverlap, K, fs, 'yaxis');
title('Spectrogram of Original Signal')
ylim([0, 60])
subplot(1,2,2)
spectrogram(ref_input, rectwin(L), noverlap, K, fs, 'yaxis');
title('Spectrogram of Synthetic Reference input')
ylim([0, 60])
figure(14)
for k = 1:9
    subplot(3,3,k)
    spectrogram(LMS(k,:), rectwin(L), noverlap, K, fs, 'yaxis');
    if k <= 3
        title(['\mu = ', num2str(mu(k)),' M = ', num2str(M(1))])
    elseif k >= 7
        title(['\mu = ', num2str(mu(k-6)),' M = ', num2str(M(3))])
    else
       title(['\mu = ', num2str(mu(k-3)),' M = ', num2str(M(2))])
    end
    ylim([0, 60])
end

figure(15)
[pxx1,f1] = periodogram(POz,[],16382,1200);
plot(f1,10*log10(pxx1))
hold on
[pxx1,f1] = periodogram(LMS(5,:),rectwin(length(POz)),16382,1200);
plot(f1,10*log10(pxx1))  
xlim([0 60])
xlabel('Frequency (Hz)')
ylabel('Periodogram')
title('Periodogram of corrupted signal and de-noised signal (M = 18 and \mu = 0.005')
legend('Corrupted signal','De-noised signal')