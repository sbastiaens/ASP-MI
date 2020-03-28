
function LMS = ALE(s, delay, M)
   LMS = zeros(100,1000);
   e = zeros(100,1000);
   mu = 0.005;
        for k=1:100
            clear u
             w = zeros(M,1);
            for n=M+delay:1000
                for i=1:M
                    u(i,n) = s(k,n - delay - i + 1);
                end
                LMS(k,n) = w'*u(:,n);
                e(k,n) = s(k,n) - LMS(k,n);
                w = w + mu*e(k,n)*u(:,n); 
            end
        end
end