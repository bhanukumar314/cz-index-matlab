
% Convert STM path to reduced STM path in symplectic basis
% Translated to MATLAB by Bhanu Kumar from original Python function written by Otto van Koert

function s_path = get_sympl_path(stm_integrator_output, mu0) %outU1, outV1, outU2, outV2)
    s_path = zeros(length(stm_integrator_output),4,4);
    pt = stm_integrator_output(1,1:6);
    [~, ~, U1_0, V1_0, U2_0, V2_0] = get_symplectic_frame(pt, mu0);

    for i = 1:length(stm_integrator_output)
        pt = stm_integrator_output(i,1:6);
        stm = reshape(stm_integrator_output(i,7:end),6,6);
        u1 = stm*U1_0; u2 = stm*U2_0; v1 = stm*V1_0; v2 = stm*V2_0; 

        [~, ~, U1, V1, U2, V2] = get_symplectic_frame(pt, mu0);
        s_path(i,:,1) = transverse_vector(u1, U1, V1, U2, V2);
        s_path(i,:,2) = transverse_vector(u2, U1, V1, U2, V2);
        s_path(i,:,3) = transverse_vector(v1, U1, V1, U2, V2);
        s_path(i,:,4) = transverse_vector(v2, U1, V1, U2, V2);
    end
end

function out_vec = transverse_vector(v, U1, V1, U2, V2)
    out_vec = [sdot6(v, V1), sdot6(v, V2), -sdot6(v, U1), -sdot6(v, U2)];
end

% Maybe this one is correct, check later:
% function out_vec = transverse_vector(v, U1, V1, U2, V2)
%     out_vec = [sdot6(V1,v), sdot6(V2,v), -sdot6(U1,v), -sdot6(U2,v)];
% end


% 
% Compute symplectic basis at a point
% 
function [Z0, XH0, U1, V1, U2, V2] = get_symplectic_frame(pt, mu0)
    I0 = zeros(6);
    I0(1:2,4:5) = -eye(2);
    I0(4:5,1:2) = eye(2);

    J0 = zeros(6);
    J0(1,2) = -1;
    J0(2,1) = 1;
    J0(4,5) = 1;
    J0(5,4) = -1;

    K0 = I0 * J0;

    XH0 = crtbp_ham(0, pt, mu0);
    DH0 = [-XH0(4), -XH0(5), -XH0(6), XH0(1), XH0(2), XH0(3)]';
    sympl_inpr = sdot6(DH0, XH0);
    Z0 = DH0 ./ sympl_inpr;
    nDH = dot(DH0, DH0);
    
    U1 = J0 * DH0;
    inpr = dot(DH0, U1);
    U1 = U1 - (inpr / nDH) * DH0;
    sympl_inpr = sdot6(Z0, U1);
    U1 = U1 - sympl_inpr * XH0;

    V1 = K0 * DH0;
    V1 = V1 - (dot(DH0, V1) / nDH) * DH0;
    V1 = V1 - sdot6(Z0, V1) * XH0;
    U1 = U1 ./ sdot6(U1, V1);
    
    U2 = [0, 0, 0, 0, 0, 1.0]';
    V2 = [0, 0, 1.0, 0, 0, 0]';

    U2 = U2 - dot(DH0, U2) / nDH * DH0;
    U2 = U2 - sdot6(Z0, U2) * XH0;
    U2 = U2 - sdot6(U1, U2) * V1;
    U2 = U2 + sdot6(V1, U2) * U1;

    V2 = V2 - dot(DH0, V2) / nDH * DH0;
    V2 = V2 -  sdot6(Z0, V2) * XH0;
    V2 = V2 -  sdot6(U1, V2) * V1;
    U2 = U2 + sdot6(V1, V2) * U1;
    U2 = U2 ./ sdot6(U2, V2);

    % frame = [Z0, XH0, U1, V1, U2, V2];
end

% 
% Variables and functions for symplectic forms and pairing
% 
function J0 = get_omega0(dim)
    J0 = zeros(dim);
    J0(dim/2+1:dim,1:dim/2) = eye(dim/2);
    J0(1:dim/2,dim/2+1:dim) = -eye(dim/2);
end

function val = sdot6(v, w)
    Omega6 = get_omega0(6);
    val = dot(v, Omega6 * w );
end
