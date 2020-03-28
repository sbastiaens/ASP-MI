
discard = 500;

N = 500;
sigma = .25;
a1 = 0.1;
a2 = 0.8;
a = [1 -a1 -a2];
mu = 0.05;
error = zeros(100,N);
w_history = zeros(2,N);


%% LMS Filter
for k = 1:100
    %setup
    w = sigma*randn(1,N+discard);
    x = filter(1,a,w)';
    x = x(discard+1:end);
    w_e = zeros(length(x),2);
    %init error vector
    e = zeros(1,N);
    e(1) = x(1);
    
    for i=2:length(x)
        e(i) = x(i) - w_e(i,:)*[x(i) x(i-1)]';
        w_e(i+1,:) = w_e(i,:)' + mu*[x(i) x(i-1)]' * e(i);
    end    
    error(k,:) = e;
end
figure,plot(1:500, 10*log10(error(1,:)))