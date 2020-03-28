
function LMS = ANC(s, epsilon, M)
   noise_est = zeros(100,1000);
   LMS = zeros(100,1000);
   mu = 0.01;
    for k=1:100
        clear u
        w = zeros(M,1); 
        for n=1:M-1  
         LMS(k,n) = s(k,n);
        end
        for n=M:1000
            for i=1:M
                u(i,n) = epsilon(k,n - i + 1);
            end
            noise_est(k,n) = w'*u(:,n);
            LMS(k,n) = s(k,n) - noise_est(k,n);
            w = w + mu*LMS(k,n)*u(:,n); 
        end
    end
end