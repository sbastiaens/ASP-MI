function [LMS,e,w_ad]=lms(x,z,mu,order)
    N=length(x); 
    w = zeros(order,1);
    w_ad = zeros(order,N);
    LMS = zeros(1,N);
    n=1;
    for i=order:N
        xsum = x(i:-1:i-order+1);
       LMS(i) = w'*xsum; 
       e(i) = z(i) - LMS(i);
       w = w + mu*e(i)*xsum; 
       w_ad(:,n) = w;
       n = n+1;
      %mu = mu - 0.0001; %for 4.3 time varying
    end 
end