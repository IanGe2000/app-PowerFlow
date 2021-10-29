%function name: NewtonRaphson
%               Newton-Raphson Method of PowerFlow Calculation

function [convergence, U_r, cita_r, Pi_r, Qi_r, s_r, S_r, U_log, cita_log] = NewtonRaphson(node, connection, Tmax, limit)
    mute = false; % mute the debug displays

    n = size(node,1);                                                       % n = number of nodes
    nPQ = 0;                                                                % nPQ = number of PQ nodes
    nPV = 0;                                                                % nPV = number of PV nodes
    nVcita = 0;                                                             % nPV = number of Vcita nodes
    
    realnameVcita = 0;                                                            % NoVcita = the name of the only Vcita node
    
    node_x = zeros(n, 6);                                                   % generates the name tables
    
    % write node_x                                                          % node_x format: [PQ nodes; PV nodes; Vcita nodes]  
    node_x(:,1) = node(:,1);                                                % node_x format: [nickname, realname, type, P1, P2, vol.lv.]
    for i = 1:n
        if node(i,2) == 1
            nPQ = nPQ + 1;
            node_x(nPQ,2:6) = node(i,1:5);
        end
    end
    for i = 1:n
        if node(i,2) == 2
            nPV = nPV + 1;
            node_x(nPQ + nPV,2:6) = node(i,1:5);
        end
    end
    for i = 1:n
        if node(i,2) == 3
            nVcita = nVcita + 1;
            realnameVcita = i;
            node_x(nPQ + nPV + nVcita,2:6) = node(i,1:5);
        end
    end
    
    % write name tables
    realnametable = node_x(:,1:2);                                          % realnametable format: [nickname(ascending), realname]
    nicknametable = zeros(n,2);                                             % nicknametable format: [nickname, realname(ascending)]
    for i = 1:n
        nicknametable(node_x(i,2),1:2) = node_x(i,1:2);
    end
    
    nicknameVcita = nicknametable(realnameVcita,1);
    
    line = zeros(1,4);                                                      % raw info container: line
    trans_T = zeros(1,4);                                                   % raw info container: transformer Type-T parameters
    trans_6 = zeros(1,10);                                                  % raw info container: transformer with 6 parameters 
    
    nline = 0;                                                              % nline = number of transmission lines
    ntrans = 0;                                                             % ntrans = number of transformers
    
    % write line, trans_T, trans_6, trans_pi
    for i = 1:size(connection,1)
        if        (connection(i,3) == 1)
            nline = nline + 1;
            line(nline,1:2) = connection(i,1:2);
            line(nline,3) = connection(i,4) + 1i * connection(i,5);
            line(nline,4) = 1i * connection(i,6);                           % line format: [i, j, R+jX, Y/2]
        else
            ntrans = ntrans + 1;
            if    (connection(i,3) == 2)
                trans_T(ntrans,1:2) = connection(i,1:2);
                trans_T(ntrans,3) = connection(i,4) + 1i * connection(i,5);
                trans_T(ntrans,4) = connection(i,6);                        % trans_T format: [i, j, R+jX, k]
            else%  connection(i,3) == 3
                trans_6(ntrans,:) = connection(i,:);
                trans_T = Trans62T(trans_6);
            end
            trans_pi = TransT2pi(trans_T);                                  % trans_pi format: [i, j, Zij, Yi0, Yj0]
        end
    end
    
    % write line_x, trans_pi_x
    if (nline ~= 0)
        line_x = zeros(nline,4);                                            % line_x format: [nickname of i, nickname of j, R+jX, k]
        for i = 1:nline
            for j = 1:2
                line_x(i,j) = nicknametable(line(i,j),1);
            end
        end
        line_x(:,3:4) = line(:,3:4);
    end
    if (ntrans ~= 0)
        trans_pi_x = zeros(ntrans,5);                                       % trans_pi_x format: [nickname of i, nickname of j, Zij, Yi0, Yj0]
        for i = 1:ntrans
            for j = 1:2
                trans_pi_x(i,j) = nicknametable(trans_pi(i,j),1);
            end
        end
        trans_pi_x(:,3:5) = trans_pi(:,3:5);
    end
    
    P = zeros(1,n);                                                         % P = active power flowing into the nodes
    Q = zeros(1,n);                                                         % Q = reactive power flowing into the nodes
    U = node_x(:,6)';                                                       % U = initial value of node voltage magnitude
    cita_deg = zeros(1,n);                                                  % cita = initial value of node voltage phase angle
    
    % write P, Q, U, cita
    for i=1:n
        if node_x(i,3) == 1
            P(node_x(i,1)) = node_x(i,4);
            Q(node_x(i,1)) = node_x(i,5);
        end
        if node_x(i,3) == 2
            P(node_x(i,1)) = node_x(i,4);
            U(node_x(i,1)) = node_x(i,5);
        end
        if node_x(i,3) == 3
            U(node_x(i,1)) = node_x(i,4);
            cita_deg(node_x(i,1)) = node_x(i,5);
        end
    end
    cita_rad = deg2rad(cita_deg);                                           % convert deg into rad
    
    U_log = U;
    cita_log = cita_deg;

    % display initial value
    if (~mute)
        disp('U = ');
        disp(U);
        disp('cita = ');
        disp(cita_deg);
    end
    
    Y = zeros(n);                                                           % Y = bus admittance matrix
    y = zeros(n);                                                           % y = admittance adjacency matrix
    
    % write y
    % write off-diagonal elements of y: admittance between nodes
    % it is assumed that the same pair of nodes cannot be connected by a
    % transmission line and a transformer at the same time
    if (nline ~= 0)
        for i = 1:size(line_x,1)
            ii = line_x(i,1); jj = line_x(i,2);
            y(ii,jj) = 1/line_x(i,3);
            y(jj,ii) = y(ii,jj);
        end
    end
    if (ntrans ~= 0)
        for i = 1:size(trans_pi_x,1)
            ii = trans_pi_x(i,1); jj = trans_pi_x(i,2);
            y(ii,jj) = 1/trans_pi_x(i,3);
            y(jj,ii) = y(ii,jj);
        end
    end
    % write diagonal elements of y: admittance between the node and ground
    if (nline ~= 0)
        for i = 1:size(line_x,1)
            ii = line_x(i,1); jj = line_x(i,2);
            y(ii,ii) = y(ii,ii) + line_x(i,4);
            y(jj,jj) = y(jj,jj) + line_x(i,4);
        end
    end
    if (ntrans ~= 0)
        for i = 1:size(trans_pi_x,1)
            ii = trans_pi_x(i,1); jj = trans_pi_x(i,2);
            y(ii,ii) = y(ii,ii) + trans_pi_x(i,4);
            y(jj,jj) = y(jj,jj) + trans_pi_x(i,5);
        end
    end
    
    % display y
    if (~mute)
        disp('y = ');
        disp(y);
    end
    
    % write Y
    for i = 1:n
        for j = 1:n
            if i == j
                Y(i,j) = sum(y(:,i));
            else
                Y(i,j) = -y(i,j);
            end
        end
    end
    
    G = real(Y);                                                            % real part of Y
    B = imag(Y);                                                            % imagine part of Y
    
    % display Y
    if (~mute)
        disp('Y = ');
        disp(Y);
    end
    
    [dP, dQ, Pi, Qi] = imbalance(n, nPQ, P, Q, U, G, B, cita_rad);          % 1st calculation of imbalance
    
    convergence = 0;                                                        % preset convergence as 0(divergent)
    for i = 1:Tmax                                                          % iteration loop
        % display iteration info
        if (~mute)
            fprintf('Iteration %d:\n', i);
        end
        
        J = Jacobi(n, nPQ, U, cita_rad, B, G, Pi, Qi);
        
        % display J
        if (~mute)
            disp('J = ');
            disp(J);
        end
        
        [dU, dcita_rad] = correction(n, nPQ, U, dP, dQ, J);
        
        % display correction
        if (~mute)
            disp('dU = ');
            disp(dU);
            disp('dcita = ');
            disp(rad2deg(dcita_rad));
        end

        U = U + dU;
        cita_rad = cita_rad + dcita_rad;
        
        U_log = [U_log; U];
        cita_log = [cita_log; rad2deg(cita_rad)];
        
        [dP, dQ, Pi, Qi] = imbalance(n, nPQ, P, Q, U, G, B, cita_rad);

        if (max(abs(dU)) < limit)
            convergence = 1;                                                % set convergence as 1(convergent)
            
            %display convergence info
            if (~mute)
                fprintf('Convergent after %d iteration(s)!\n', i);
            end%if
            
            break%for;
        end%if
        
    end%for
        
    % Calculate power flow
    if (convergence)
        s = zeros(n);                                                       % the power lost on each connection (directional)
        S = zeros(n);                                                       % the power lost on each connection (sum of two direction)
        S_loss_1 = 0;                                                       % the total power loss (calculation method 1)
        S_loss_2 = 0;                                                       % the total power loss (calculation method 2)
        
        % write s
        if (nline ~= 0)
            for i = 1:size(line_x,1)
                ii = line_x(i,1); jj = line_x(i,2);
                Uii = U(ii)*exp(1i*cita_rad(ii)); Uiiconj = conj(Uii);
                Ujj = U(jj)*exp(1i*cita_rad(jj)); Ujjconj = conj(Ujj);
                s(ii,jj) = Uii * (Uiiconj * conj(line_x(i,4)) + (Uiiconj - Ujjconj)*conj(y(ii,jj)));
                s(jj,ii) = Ujj * (Ujjconj * conj(line_x(i,4)) + (Ujjconj - Uiiconj)*conj(y(jj,ii)));
            end
        end
        if (ntrans ~= 0)
            for i = 1:size(trans_pi_x,1)
                ii = trans_pi_x(i,1); jj = trans_pi_x(i,2);
                Uii = U(ii)*exp(1i*cita_rad(ii)); Uiiconj = conj(Uii);
                Ujj = U(jj)*exp(1i*cita_rad(jj)); Ujjconj = conj(Ujj);
                s(ii,jj) = Uii * (Uiiconj * conj(trans_pi_x(i,4)) + (Uiiconj - Ujjconj)*conj(y(ii,jj)));
                s(jj,ii) = Ujj * (Ujjconj * conj(trans_pi_x(i,5)) + (Ujjconj - Uiiconj)*conj(y(jj,ii)));
            end
        end
        
        % write S (adding the symmetrical elements of s to construct a upper triangular matrix)
        for i = 1:n
            for j = 1:n
                if (i < j)
                    S(i,j) = S(i,j) + s(i,j);
                elseif (i > j)
                    S(j,i) = S(j,i) + s(i,j);
                else
                    S(i,j) = s(i,j);
                end
            end
        end
        
        % write S_loss_1
        for i = 1:n
            if (i ~= nicknameVcita)
                S_loss_1 = S_loss_1 + node_x(i,4) + 1i*node_x(i,5);
            end
        end
        S_loss_1 = S_loss_1 + U(nicknameVcita)*exp(1i*cita_rad(nicknameVcita)) * conj(U.*exp(1i*cita_rad)) * conj(Y(:,nicknameVcita));
        
        % write S_loss_2
        S_loss_2 = sum(sum(s));
        
        % display s, S and the sum of power consumption calculation
        if (~mute)
            fprintf("s = \n");
            disp(s);
            fprintf("S = \n");
            disp(S);
            fprintf("total power loss 1 = ");
            disp(S_loss_1)
            fprintf("total power loss 2 = ");
            disp(S_loss_2);
            fprintf("difference between the two calcualting method = ");
            disp(S_loss_1 - S_loss_2);
        end
    else
        s = NaN;
        S = NaN;
    end
    
    U_r = zeros(1,n);                                                       % output
    cita_r = zeros(1,n);
    Pi_r = zeros(1,n);
    Qi_r = zeros(1,n);
    if (convergence)
        s_r = zeros(n);
        S_r = zeros(n);
    else
        s_r = NaN;
        S_r = NaN;
    end
    
    for i = 1:n
        U_r(realnametable(i,2)) = U(i);
        cita_r(realnametable(i,2)) = rad2deg(cita_rad(i));
        Pi_r(realnametable(i,2)) = Pi(i);
        Qi_r(realnametable(i,2)) = Qi(i);
    end%for
    
    if (convergence)
        for i = 1:n
            for j = 1:n
                s_r(realnametable(i,2),realnametable(j,2)) = s(i,j);
                S_r(realnametable(i,2),realnametable(j,2)) = S(i,j);
            end
        end%for
    end%if

end%function
