% function [cz_idx, cz_planar, cz_spatial, s_path, eigs, phi] = get_split_cz_index(state, period, mu)
%
% MATLAB function to compute CRTBP planar periodic orbit Conley-Zehnder indices.
% Helper functions all translated to MATLAB by Bhanu Kumar, from the original  
% Python library "cz-index" written by Otto van Koert (Github: ovkoert)
%
% Recommended citations (please cite both):
%
% For this MATLAB code:
% B. Kumar and A. Moreno, “Networks of Periodic Orbits in the Earth-Moon System Through a Regularized
% and Symplectic Lens,” AAS/AIAA Astrodynamics Specialist Conference, Aug 2025. Paper AAS-25-677. 
% 
% For the methodology and original Python library by Otto van Koert: 
% A. Moreno, C. Aydin, O. v. Koert, U. Frauenfelder, and D. Koh, “Bifurcation Graphs for the CR3BP via
% Symplectic Methods,” J. Astronaut. Sci., Vol. 71, No. 6, 2024, p. 51.
%
% inputs:
%  state = state for periodic orbit (6x1 or 1x6, should be a planar periodic orbit)
%  period = period of periodic orbit
%  mu = mass ratio
%
% outputs:
%  cz_idx: Conley-Zehnder index of periodic orbit
%  cz_planar: planar Conley-Zehnder index of periodic orbit
%  cz_spatial: spatial Conley-Zehnder index of periodic orbit
%  s_path: path of reduced state transition matrices used for CZ index calculation
%  eigs: eigenvalues of monodromy
%  phi: monodromy matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cz_idx, cz_planar, cz_spatial, s_path, eigs, phi] = get_split_cz_index(state, period, mu, ~, ~)

    if length(state) ~= 6 
      fprintf('error:  state must be of length 6\n');
      cz_idx = nan; cz_planar = nan; cz_spatial = nan; s_path = nan; eigs = nan*state; phi = nan*state;
      return;
    end

    if state(3) ~= 0 || state(6) ~= 0 
      fprintf('error:  periodic orbit must be planar with z and zdot = 0\n');
      cz_idx = nan; cz_planar = nan; cz_spatial = nan; s_path = nan; eigs = nan*state; phi = nan*state;
      return;
    end
    
    x0 = state;
    x0(4) = x0(4) - x0(2); %Convert to 
    x0(5) = x0(5) + x0(1); %momenta from velocity
    
    phi=eye(6,6);
    x0(7:6*6+6)=reshape(phi,6*6,1);
    options = odeset('RelTol', 2.23e-14, 'AbsTol', ones(length(x0),1)*2.23e-14);
    [~, xout] = ode45(@crtbpstm_ham,[0,period],x0,options, mu);
    
%     norm(xout(end,1:6) - xout(1,1:6))

    [cz_planar, cz_spatial, s_path] = get_split_index(xout, mu);
    cz_idx = cz_planar + cz_spatial;
    
    % Also output eigenvalues of monodromy matrix
    xfinal = xout(end,:);
    phi = reshape(xfinal(7:6*6+6),6,6);
    phi2 = phi([1 2 4 5 3 6], :); %separate planar and spatial
    phi2 = phi2(:, [1 2 4 5 3 6]); %parts in block form
    
    if all(all(phi2(5:6, 1:4) == 0)) && all(all(phi2(1:4, 5:6) == 0))
        eigs_p = eig(phi2(1:4, 1:4));
        eigs_s = eig(phi2(5:6, 5:6));
        [~, idx_p] = sort(abs(eigs_p-1.0));
        [~, idx_s] = sort(abs(eigs_s));
        eigs = [eigs_p(idx_p); eigs_s(idx_s)];
    else
        eigs = eig(phi);
    end
    
end

function [planar_CZ, out_plane_CZ, s_path] = get_split_index(stm_integrator_output, mu0)
    s_path = get_sympl_path(stm_integrator_output, mu0);
    T = size(s_path,1);
    non_split = 0.0;
    for i = 1:2
        for j =1:2
            non_split = non_split + sum(abs(s_path(:, 2*i -1, 2*j))) + sum(abs(s_path(:, 2*i, 2*j -1)) ) ;
        end
    end
    if non_split ~= 0.0
        fprintf("This symplectic path does not split: %.10f", sum(non_split))
        planar_CZ = nan; out_plane_CZ = nan; 
        return
    end
    p_path = zeros(size(s_path,1),2,2);
    p_path(:,1,1) = s_path(:,1,1);
    p_path(:,1,2) = s_path(:,1,3);
    p_path(:,2,1) = s_path(:,3,1);
    p_path(:,2,2) = s_path(:,3,3);
    planar_CZ = planar_index(p_path, T);
    o_path = zeros(size(s_path,1),2,2);
    o_path(:,1,1) = s_path(:,2,2);
    o_path(:,1,2) = s_path(:,2,4);
    o_path(:,2,1) = s_path(:,4,2);
    o_path(:,2,2) = s_path(:,4,4);
    out_plane_CZ = planar_index(o_path, T);
end

function cdets = retract_sympl_path(path)
    % Retracts a symplectic path onto U(n), and takes its determinant
    % 
    % Parameters
    % ----------
    % path : array
    % Array with shape (N,2n,2n)
    % 
    % Returns
    % -------
    % array
    % Array containing the path of complex determinants 
    % 
    N = size(path,1);
    dim = size(path,2);
    e_vals = zeros(size(path));
    for i=1:N
        [V, ~, W] = svd(reshape(path(i,:,:),dim,dim));
        e_vals(i,:,:) = V * W';
    end
    cdets = zeros(N,1);
    for i = 1:size(e_vals,1)
        cdets(i) = det(reshape(e_vals(i,1:dim/2,1:dim/2) + 1j * e_vals(i,dim/2+1:dim,1:dim/2),dim/2,dim/2) );
    end 
end
    
% Functions for dealing with paths in Sp(2)

function parity = parity_type(s_path_xi, T)
    traceVal = trace(reshape(s_path_xi(T,:,:),2,2));
    if traceVal < 2
        parity = 1;
    else
        parity = 0;
    end
end

function val = return_parity_corrected(v, parity)
    n = floor(v);
    if mod(n,2) == parity
        val = n;
    else
        val = n+1;
    end
end

function frac_idx = compute_fractional_angular_index(cdets, T)
    angle = 0;
    for i = 1:T
        delta = cdets(i+1) / cdets(i);
        angle = angle + imag(log(delta));
    end
    frac_idx = (angle / pi);
end

function planar_CZ = planar_index(p_path, T)
    cdets = retract_sympl_path(p_path);
    angular_idx = compute_fractional_angular_index(cdets, length(cdets)-1);
    fprintf("\nSum of angle jumps divided by pi (need not be near integer): %.10f \n", angular_idx);
    planar_CZ = return_parity_corrected(angular_idx, parity_type(p_path, T));
end


