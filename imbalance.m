%Function name: imbalance
%               calculate imbalance for Newton-Laphson Method

function [dP, dQ, Pi, Qi] = imbalance(n, m, P, Q, U, G, B, cita)

    Pn = zeros(1, n);                                                       %preallocation of Pn
    Pi = zeros(1, n);                                                       %preallocation of Pi
    
    for i=1:n
        for j=1:n
            Pn(j)=U(i)*U(j)*(G(i,j)*cos(cita(i)-cita(j))+B(i,j)*sin(cita(i)-cita(j)));
        end
        Pi(i)=sum(Pn);                                                      %>>Pi
    end
    
    dP=P(1:n-1)-Pi(1:n-1);                                                  %>>dP

    
    Qn = zeros(1, n);                                                       %preallocation of Qn
    Qi = zeros(1, n);                                                       %preallocation of Qi
    
    for i=1:n
        for j=1:n
            Qn(j)=U(i)*U(j)*(G(i,j)*sin(cita(i)-cita(j))-B(i,j)*cos(cita(i)-cita(j)));
        end
        Qi(i)=sum(Qn);                                                      %>>Qi
    end
    
    dQ=Q(1:m)-Qi(1:m);                                                      %>>dQ
    
end%function