function [Hstar,obj] = mySolvingHIter(Kmatrix,Hstar,H0,lambda)

flag =1;
iter = 0;
while flag
    iter = iter+1;
    M = 2*(Kmatrix*Hstar) + lambda*H0;
    [U,S,V] = svd(M,'econ');
    Hstar = U*V';
    obj(iter) = trace(Hstar'*(Kmatrix*Hstar + lambda*H0) );
  
    if(iter>2&& ((obj(iter)-obj(iter-1))/obj(iter))<1e-4)
        flag = 0;
    end
end