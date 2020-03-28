%% ASP&MI Part 4 

clear all
close all

%% 1: LMS
load('time-series.mat')
y = y - mean(y);
y = y';

mu = 1*10^(-5);
M = 4;
N = length(y);

w = zeros(M,N);
e = zeros(1,N);
yhat = zeros(1,N);
xsum = zeros(M,N);
for n = M+1:N
    xsum(:,n) = y(1,(n-1):-1:n-M);
    yhat(n) = w(:,n)'*xsum(:,n);
    e(n) = y(n) - yhat(n);
    w(:,n+1) = w(:,n) + mu*e(n)*xsum(:,n); 
end

figure(1)
plot(1:N, y)
hold on 
plot(1:N, yhat)
xlabel('Sample Number')
ylabel('Amplitude')
legend('Original', 'Estimate')

MSE = mean(e.^2);
var_output = var(yhat);
var_error = var(e);
Gain = 10*log10(var_output/var_error);

%% 2: Dynamical perceptron with tanh
w = zeros(M,N);
e = zeros(1,N);
yhat = zeros(1,N);
xsum = zeros(M,N);

for n = M+1:N
    xsum(:,n) = y(1,(n-1):-1:n-M);
    yhat(n) = tanh(w(:,n)'*xsum(:,n));
    e(n) = y(n) - yhat(n);
    w(:,n+1) = w(:,n) + mu*e(n)*xsum(:,n); 
end

figure(2)
plot(1:N, y)
hold on 
plot(1:N, yhat)
xlabel('Sample Number')
ylabel('Amplitude')
legend('Original', 'Estimate')

%% 3: Dynamical perceptron with a*tanh (a = 50)

w = zeros(M,N);
e = zeros(1,N);
yhat = zeros(1,N);
xsum = zeros(M,N);

for n = M+1:N
    xsum(:,n) = y(1,(n-1):-1:n-M);
    yhat(n) = 50*tanh(w(:,n)'*xsum(:,n));
    e(n) = y(n) - yhat(n);
    w(:,n+1) = w(:,n) + mu*e(n)*xsum(:,n); 
end

figure(3)
plot(1:N, y)
hold on 
plot(1:N, yhat)
xlabel('Sample Number')
ylabel('Amplitude')
legend('Original', 'Estimate')

MSE_1 = mean(e.^2);
var_output = var(yhat);
var_error = var(e);
Gain_1 = 10*log10(var_output/var_error);

%% 4: Bias added to model

load('time-series.mat');
y = y';

w =  zeros(M+1,N);
e = zeros(1,N);
yhat = zeros(1,N);
xsum = zeros(M+1,N);

for n = M+1:N
    xsum(:,n) = [1, y(1,(n-1):-1:n-M)];
    yhat(n) = 48*tanh((w(:,n)'*xsum(:,n)));
    e(n) = y(n) - yhat(n);
    w(:,n+1) = w(:,n) + mu*e(n)*xsum(:,n); 
end

figure(4)
plot(1:N, y)
hold on 
plot(1:N, yhat)
xlabel('Sample Number')
ylabel('Amplitude')
legend('Original', 'Estimate')

MSE_2 = mean(e.^2);
var_output = var(yhat);
var_error = var(e);
Gain_2 = 10*log10(var_output/var_error);

%% 5: Pre-trained weights method
y_20 = y(1,1:20);
N = length(y_20);
M = 4;

w = zeros(M+1,1);
e = zeros(1,N);
yest = zeros(1,N);
xsum = zeros(M+1,N);
res = zeros(M+1,100);

for i=1:100
      for n = M+1:N
        xsum(:,n) = [1, y(1,(n-1):-1:n-M)];
        yest(n) = 50*tanh((w(2:5,1)'*xsum(2:5,n))+(w(1,1)*xsum(1,n)));
        e(n) = y(n) - yest(n);
        w = w + mu*e(n)*xsum(:,n); 
      end
    res(:,i) = w;
end


N = length(y);
w_init = res(:,end);
yhat = zeros(1,N);
xsum = zeros(M+1,N);
e = zeros(1,N);

for n = M+1:N
    xsum(:,n) = [1, y(1,(n-1):-1:n-M)];
    yhat(n) = 48.78*tanh((w_init(2:5,1)'*xsum(2:5,n))+(w_init(1,1)*xsum(1,n)));
    e(n) = y(n) - yhat(n);
    w_init = w_init + mu*e(n)*xsum(:,n); 
end

figure(5)
plot(1:N, y)
hold on 
plot(1:N, yhat)
xlabel('Sample Number')
ylabel('Amplitude')
legend('Original', 'Estimate')

MSE_3 = mean(e.^2);
var_output = var(yhat);
var_error = var(e);
Gain_3 = 10*log10(var_output/var_error);

