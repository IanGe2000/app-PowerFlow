%function name: Jacobi
%               calculate the Jacobi matrix for Newton-Laphson Method

function [J] = Jacobi(n, m, U, cita, B, G, Pi, Qi)

%Target object: J
%J format:      J = [H N; K L]      (blocked)

    H = zeros(n-1, n-1);                                                    %preallocation of H
    N = zeros(n-1, m);                                                      %                 N
    K = zeros(m, n-1);                                                      %                 K
    L = zeros(m, m);                                                        %                 L

    %Calculate Off-diagonal elements
    %       i!=j
    %       diagonal elements will be recalculated
    
    for i=1:n-1
        for j=1:n-1
            H(i,j)=U(i)*U(j)*(G(i,j)*sin(cita(i)-cita(j))-B(i,j)*cos(cita(i)-cita(j)));
        end
    end
    for i=1:n-1
        for j=1:m
            N(i,j)=U(i)*U(j)*(G(i,j)*cos(cita(i)-cita(j))+B(i,j)*sin(cita(i)-cita(j)));
        end
    end
    for i=1:m
        for j=1:n-1
            K(i,j)=-U(i)*U(j)*(G(i,j)*cos(cita(i)-cita(j))+B(i,j)*sin(cita(i)-cita(j)));
        end
    end
    for i=1:m
        for j=1:m
            L(i,j)=U(i)*U(j)*(G(i,j)*sin(cita(i)-cita(j))-B(i,j)*cos(cita(i)-cita(j)));
        end
    end
    
    %Re-Calculate Diagonal elements
    %       i=j
    
    for i=1:n-1
        H(i,i)=-U(i).^2*B(i,i)-Qi(i);
    end
    for i=1:m
        N(i,i)=U(i).^2*G(i,i)+Pi(i);
    end
    for i=1:m
        K(i,i)=-U(i).^2*G(i,i)+Pi(i);
    end
    for i=1:m
        L(i,i)=-U(i).^2*B(i,i)+Qi(i);
    end

    J=[H N; K L];                                                           %>>J

end%function
