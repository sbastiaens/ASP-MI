function [yhat, h, e]= DTFCLMS(x,y,M,N,mu,gamma)
    h = zeros(M,N);
    e = zeros(1,N);
    yhat = zeros(1,N);
        for n = 1:N
            yhat(n) = h(:,n)'*x(:,n);
            e(n) = y(n) - yhat(n);
            h(:,n+1) =(1- mu*gamma)*h(:,n) + mu*conj(e(n))*x(:,n);
        end
        h = h(:, 2:end);
end