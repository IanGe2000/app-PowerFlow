%function name: Get_Trans
%               user input manipulation

function [Trans_T] = Trans62T(Trans_6)

%Target object: Trans_6
%Trans_6 format: [i, j, 3, Sn, U1n, U2n, Pk, Uk(%), Po, Io(%)]
%Trans format: [i, j, R+jX, k]
%               with Zt=R+jX being the impedance reduced to node j side

% IMPORTANT: This convertion results in the following model, notice that
% the impedace Zt is reduced to the node j side
%
%                       i    k:1            j
%                       |----OO----[ Zt ]---|
%                       |                   |
%


Trans_T = zeros(size(Trans_6,1),4);

for a = 1:size(Trans_6,1)
    Trans_T(a,1) = Trans_6(a,1);                                              % Write Node i
    Trans_T(a,2) = Trans_6(a,2);                                              % Write Node j
    Rt = Trans_6(a,7) / (1000*Trans_6(a,4)) * (Trans_6(a,6)^2 / Trans_6(a,4));% Calculate Rt
    Xt = Trans_6(a,8) / (100) * (Trans_6(a,6)^2 / Trans_6(a,4));              % Calculate Xt
    Trans_T(a,3) = Rt + 1i * Xt;                                              % Write R+iX
    Trans_T(a,4) = (Trans_6(a,5) / Trans_6(a,6));                             % Write k
    
end%for

end%function
    