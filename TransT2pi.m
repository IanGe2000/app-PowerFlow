%function name: Get_Trans_pi
%               user input manipulation

function [Trans_pi] = TransT2pi(Trans_T)

% Target object: Trans
% Trans format: [i, j, Zt(=R+jX), k]
% Trans_pi format: [i, j, Zij, Yi0, Yj0]

% IMPORTANT: This convertion takes the following model as an input, notice
% that the impedace Zt is reduced to the node j side
%
%                       i    k:1            j
%                       |----OO----[ Zt ]---|
%                       |                   |
%
% and the result of this convertion gives the impedance between the nodes
% and the admittance of between each node and the ground
%
%                       i                   j
%                       |--+---[ Zij ]---+--|
%                       |  |             |  |
%                        [Yi0]         [Yj0] 
%                          |             |
%                          G             G
%

    Zt_pi = zeros(size(Trans_T,1),3);
    Zt_pi(:,1) = Trans_T(:,3).*Trans_T(:,4);
    Zt_pi(:,2) = Trans_T(:,3).*(Trans_T(:,4)).^2./(1-Trans_T(:,4));
    Zt_pi(:,3) = Trans_T(:,3).*(Trans_T(:,4))./(Trans_T(:,4)-1);
    
    Trans_pi = [Trans_T(:,1:2) Zt_pi(:,1) 1./Zt_pi(:,2) 1./Zt_pi(:,3)];


end%function