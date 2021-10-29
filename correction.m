%function name: correction
%               calculate the correction for Newton-Laphson Method

function [dU, dcita] = correction(n, m, U, dP, dQ, J)

    Ud2 = zeros(m);                                                         %preallocation of Ud2
    
    for i=1:m
        Ud2(i,i)=U(i);
    end
    
    dPQ=[dP dQ]';
    
    dUcita=(J\dPQ)';
    
    dcita=dUcita(1:n-1);
    dcita=[dcita 0];                                                        %>>dcita
    
    dU=(Ud2*dUcita(n:n+m-1)')';
    dU=[dU zeros(1,n-m)];                                                   %>>dU

end
 