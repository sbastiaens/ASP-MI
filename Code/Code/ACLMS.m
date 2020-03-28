function [yhat, h, g, e]= ACLMS(x,y,M,N,mu,mode)
    h = zeros(M,N);
    g = zeros(M,N);
    e = zeros(1,N);
    yhat = zeros(1,N);
    xsum = zeros(M,N);
    if (strcmp(mode,'Identification'))
        for n = M:N
            xsum(:,n) = [x(n); x(n-1)];
            yhat(n) = h(:,n)'*xsum(:,n) + g(:,n)'*conj(xsum(:,n));
            e(n) = y(n) - yhat(n);
            h(:,n+1) = h(:,n) + mu*conj(e(n))*xsum(:,n);
            g(:,n+1) = g(:,n) + mu*conj(e(n))*conj(xsum(:,n));
        end
    elseif (strcmp(mode,'Prediction'))
        for n = M+1:N
            xsum(:,n) = x(1,(n-1):-1:n-M);
            yhat(n) = h(:,n)'*xsum(:,n) + g(:,n)'*conj(xsum(:,n));
            e(n) = y(n) - yhat(n);
            h(:,n+1) = h(:,n) + mu*conj(e(n))*xsum(:,n);
            g(:,n+1) = g(:,n) + mu*conj(e(n))*conj(xsum(:,n));
        end
    end
end