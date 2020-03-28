function term = gassTerm(i,k,xsum,mu,e,oldterm,method)
     if (strcmp(method,'Benveniste'))
         term = (eye(1)- mu(i-1)*xsum(2)'*xsum(2))*oldterm+e(k,i-1)*xsum(2)';
     elseif (strcmp(method,'Ang'))
         term = 0.8*oldterm+e(k,i-1)*xsum(2);
     elseif (strcmp(method,'Matthews'))
         term = e(k,i-1)*xsum(2);
     end
 end

 