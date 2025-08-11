% function [cz_idx, s_path, extend_s_path, eigs, phi] = get_cz_index(state, period, mu, steps, error_report)
%
% MATLAB function to compute CRTBP periodic orbit Conley-Zehnder indices.
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
%  state = state for periodic orbit (6x1 or 1x6)
%  period = period of periodic orbit
%  mu = mass ratio
%  steps = number of discrete steps used to construct extension path from monodromy matrix to base matrix
%  error_report = verbose reporting of errors in various quantities (if true)
%
% outputs:
%  cz_idx: Conley-Zehnder index of periodic orbit
%  s_path: path of reduced state transition matrices used for CZ index calculation
%  extend_s_path: path connecting final reduced monodromy matrix to appropriate base symplectic matrix 
%  eigs: eigenvalues of monodromy
%  phi: monodromy matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cz_idx, s_path, extend_s_path, eigs, phi] = get_cz_index(state, period, mu, steps, error_report)

    if length(state) ~= 6
      fprintf('error:  state must be of length 6\n');
      cz_idx = nan; eigs = nan*state; phi = nan*state;
      return;
    end

    if nargin < 4
      steps = 20000;
      error_report = false;
    elseif nargin < 5
      error_report = false;
    end
    
    x0 = state;
    x0(4) = x0(4) - x0(2); %Convert to 
    x0(5) = x0(5) + x0(1); %momenta from velocity
    
    phi=eye(6,6);
    x0(7:6*6+6)=reshape(phi,6*6,1);
    options = odeset('RelTol', 2.23e-14, 'AbsTol', ones(length(x0),1)*2.23e-14);
    [~, xout] = ode45(@crtbpstm_ham,[0,period],x0,options, mu);
    
%     norm(xout(end,1:6) - xout(1,1:6))

    s_path = get_sympl_path(xout, mu);
    [cz_idx, extend_s_path] = get_index_sympl_path(s_path, steps, error_report);
    
    % Also output eigenvalues of monodromy matrix
    xfinal = xout(end,:);
    phi = reshape(xfinal(7:6*6+6),6,6);
    phi2 = phi([1 2 4 5 3 6], :); %separate planar and spatial
    phi2 = phi2(:, [1 2 4 5 3 6]); %parts in block form (if planar orbit)
    
    if all(all(phi2(5:6, 1:4) == 0)) && all(all(phi2(1:4, 5:6) == 0)) % if planar orbit
        eigs_p = eig(phi2(1:4, 1:4));
        eigs_s = eig(phi2(5:6, 5:6));
        [~, idx_p] = sort(abs(eigs_p-1.0));
        [~, idx_s] = sort(abs(eigs_s));
        eigs = [eigs_p(idx_p); eigs_s(idx_s)];
    else
        eigs = eig(phi);
    end
    
end

% 
% Variables and functions for symplectic forms and pairing
% 
function J0 = get_omega0(dim)
    J0 = zeros(dim);
    J0(dim/2+1:dim,1:dim/2) = eye(dim/2);
    J0(1:dim/2,dim/2+1:dim) = -eye(dim/2);
end

function val = sdot4(v, w)
    Omega4 = get_omega0(4);
    val = dot(v, Omega4 * w );
end

function val = sdot6(v, w)
    Omega6 = get_omega0(6);
    val = dot(v, Omega6 * w );
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


function [K, A, N] = iwasawa(S)
    % Returns Iwasawa decomposition of a symplectic matrix, following 
    % 
    % Parameters
    % ----------
    % S : array
    % Array with shape (2n,2n)
    % 
    % Returns
    % -------
    % array
    % array with shape (2n,2n) containing K-part
    % array
    % array with shape (2n,2n) containing A-part
    % array
    % array with shape (2n,2n) containing N-part
    % 
    dim = size(S,1);
    n = dim / 2;
    S1 = S(:, 1:n);
    [Q, R] = qr(S1,0);
    H = zeros(n);
    for i=1:n
        H(i, i) = 1 / R(i, i);
    end
    U = H * R;
    D = zeros(n);
    sD = zeros(n);
    sD_inv = zeros(n);
    for i=1:n
        H(i, i) = R(i, i);
        D(i, i) = R(i, i) * R(i, i);
        sD(i,i) = sqrt(D(i, i));
        sD_inv(i,i) = 1 / sD(i,i);
    end
    A = zeros(dim);
    A(1:n,1:n) = sD ;
    A(n+1:dim,n+1:dim) = sD_inv;
    K = zeros(dim);
    K(1:n, 1:n) = Q(1:n,1:n) * H * sD_inv ;
    K(1:n, n+1:dim) = -Q(n+1:dim,1:n) * H * sD_inv;
    K(n+1:dim, n+1:dim) = K(1:n, 1:n);
    K(n+1:dim, 1:n) = -K(1:n, n+1:dim);
    N = zeros(dim);
    N(1:n, 1:n) = U ;
    N(:, n+1:dim) = A \ (K' * S(:,n+1:dim));
end


function U = KtoU(K)
    % Converts 2n by 2n U(n) matrix to n by n complex form
    % 
    % Parameters
    % ----------
    % path : array
    % Array with shape (2n,2n)
    % 
    % Returns
    % -------
    % array
    % array with shape (n,n) containing U(n) matrix
    % 
    n = length(K) / 2;

    U = K(1:n,1:n) - 1j * K(n+1:2*n,1:n);
end


function K = UtoK(U)
    % Converts n by n U(n) matrix to 2n by 2n real form
    % 
    % Parameters
    % ----------
    % path : array
    % Array with shape (n,n)
    % 
    % Returns
    % -------
    % array
    % array with shape (2n,2n) containing U(n) matrix
    % in real form
    % 
    n = length(U);
    K = zeros(2*n);
    K(1:n,1:n) = real(U);
    K(n+1:2*n,n+1:2*n) = real(U);
    K(n+1:2*n,1:n) = -imag(U);
    K(1:n,n+1:2*n) = +imag(U);
end

function maslov_val = maslov_func(S)
    % Returns the value of the Maslov function det(id-S)
    % 
    % Parameters
    % ----------
    % S : array
    % Array with shape (2n,2n)
    % 
    % Returns
    % -------
    % float
    % 
    maslov_val = det(eye(4)-S);
end




% Functions for dealing with elliptic matrices
% These are symplectic matrices with eigenvalues on the unit circle.
% 
% 
%  n by n U(n) matrix to 2n by 2n real form
% 
% Parameters
% ----------
% path : array
%     Array with shape (n,n)
% 
% Returns
% -------
% array
%     array with shape (2n,2n) containing U(n) matrix
%     in real form
% 

function [c_pairs, eigvals, B] = find_elliptic_pairs(S)
    % Returns the elliptic pairs of eigenvalues of a symplectic matrix S.
    % Since elliptic pairs are supposed to be on the unit circle and not degenerate,
    % we discard values that are further away than eps_tol
    % 
    % Parameters
    % ----------
    % S : array
    % Array with shape (2n,2n)
    % 
    % Returns
    % -------
    % Tuple consisting of a list with the indices of the elliptic pairs,
    % an array with the eigenvalues, and 
    % an array with the basis change (from diagonalization)
    % 
    eps_tol = 1e-5;

    [B, eigvals] = eig(S);
    eigvals = diag(eigvals);
    c_pairs = [];
    indices = [];
    for i=1:length(eigvals)
        eigval = eigvals(i);
        if abs(eigval-1) < eps_tol || abs(eigval+1)< eps_tol || abs(abs(eigval)-1) > eps_tol || any(indices == i)
            continue
        end
        [~, j] = min( abs( eigval-conj(eigvals(i+1:end)) ) );
        inv_idx = i + j;
        indices = [indices i inv_idx];
        c_pairs = [c_pairs [i;inv_idx] ];
        % for j = (i+1):length(eigvals)
        %     if abs( conj(eigvals(j)) - eigval ) < eps_tol
        %         indices = [indices i j];
        %         c_pairs = [c_pairs [i;j] ];
        %         break
        %     end
        % end
    end
end

function path_out = path_removing_elliptic_pairs(S, steps)
    % Returns path of symplectic matrices starting at S and ending at S', where
    % S' has no elliptic pairs left
    % 
    % Parameters
    % ----------
    % S : array
    % Array with shape (2n,2n)
    % steps : int (optional argument)
    % 
    % Returns
    % -------
    % path : array with shape (steps, 2n, 2n)
    % 
    
    if nargin < 2
      steps = 200;
    end
    
    [c_pairs, eigvals, B] = find_elliptic_pairs(S);
    path = zeros(steps, 4, 4);
    if isempty(c_pairs)
        for i = 1:steps
            path(i,:,:) = S;
        end
        path_out = path;
        return
    end
    
    Binv = inv(B);
    angle_paths = [];
    for i=1:size(c_pairs,2)
        pair = c_pairs(:,i);
        eigval = eigvals(pair(1));
        angle = atan2(imag(eigval), real(eigval));
        angle_paths = [angle_paths; exp(linspace(angle * 1j, pi * 1j, steps))];
    end
    
    for s = 1:steps
        lst = eigvals;
        for idx = 1:size(c_pairs,2)
            pair = c_pairs(:,idx);
            lst(pair(1)) = angle_paths(idx,s);
            lst(pair(2)) = conj(angle_paths(idx,s));
        end
        diag_s = diag(lst);
        path(s,:,:) = B * diag_s * Binv;
    end
    path_out = real(path);
end

% Functions for dealing with hyperbolic matrices
% These are symplectic matrices with eigenvalues in norm greater or less than 1,
% 

function [c_pairs, eigvals, B] = find_real_hyperbolic_pairs(S)
    % Returns the pairs of real eigenvalues of a symplectic matrix S.
    % Eigenvalues with imaginary part larger than eps, or eigenvalues that are closer
    % than eps to 1 are discarded.
    % 
    % Parameters
    % ----------
    % S : array
    % Array with shape (2n,2n)
    % 
    % Returns
    % -------
    % Tuple consisting of a list with the indices of real hyperbolic pairs,
    % an array with the eigenvalues, and 
    % an array with the basis change (from diagonalization)
    % 
    eps_tol = 1e-5;

    [B, eigvals] = eig(S);
    eigvals = diag(eigvals);

    neg_one_idxs = find(abs(eigvals+1)<eps_tol);
    % In case of eigenvalues -1, -1 creating nearly dependent eigenvectors 
    if length(neg_one_idxs) == 2
        idx_1 = neg_one_idxs(1);
        idx_2 = neg_one_idxs(2);
        % If -1 eigenvals/vectors are complex conjugates w/tiny imag part
        if norm(conj(B(:,idx_1)) - B(:,idx_2)) < eps_tol || norm(conj(B(:,idx_1)) + B(:,idx_2)) < eps_tol
            eigvec = B(:,idx_1);
            B(:,idx_1) = real(eigvec)./norm(real(eigvec));
            B(:,idx_2) = imag(eigvec)./norm(imag(eigvec));
        % If -1 eigenvectors are real but nearly dependent
        elseif norm(B(:,idx_1) - B(:,idx_2)) < 1e-1
            B(:,idx_2) = (B(:,idx_1) - B(:,idx_2)) ./ norm(B(:,idx_1) - B(:,idx_2));
        elseif norm(B(:,idx_1) + B(:,idx_2)) < 1e-1
            B(:,idx_2) = (B(:,idx_1) + B(:,idx_2)) ./ norm(B(:,idx_1) + B(:,idx_2));
        end
    end
    
    real_idxs = abs(imag(eigvals))<eps_tol;
    eigvals(real_idxs) = real(eigvals(real_idxs));
    c_pairs = [];
    indices = [];
    for i = 1:length(eigvals)
        eigval = eigvals(i);
        if abs(real(eigval) - eigval) > eps_tol || abs(eigval+1)< eps_tol || abs(eigval-1)< eps_tol || any(indices == i)
            continue
        end

        [~, j] = min( abs(eigval-1./eigvals(i+1:end)) );
        inv_idx = i + j;
        indices = [indices i inv_idx];
        c_pairs = [c_pairs [i;inv_idx] ];
    end
end

function path_out = path_removing_negative_hyperbolic_pairs(S, steps)
    % Returns path of symplectic matrices starting at S and ending at S', where
    % S' has no negative real hyperbolic pairs left
    % 
    % Parameters
    % ----------
    % S : array
    % Array with shape (2n,2n)
    % steps : int
    % 
    % Returns
    % -------
    % path : array with shape (steps, 2n, 2n)
    % 

    if nargin < 2
      steps = 200;
    end
    
    [c_pairs, eigvals, B] = find_real_hyperbolic_pairs(S);
    % eigvals(abs(imag(eigvals))<eps_tol) = real(eigvals(abs(imag(eigvals))<eps_tol));
    path = zeros(steps, 4, 4);
    if isempty(c_pairs)
        for i = 1:steps
            path(i,:,:) = S;
        end
        path_out = path;
        return
    end
    
    eigenval_paths = [];
    Binv = inv(B);
    for i=1:size(c_pairs,2)
        pair = c_pairs(:,i);
        eigval = eigvals(pair(1));
        if eigval < 0
            eigenval_paths = [eigenval_paths; linspace(eigval, -1.0, steps)];
        else
            eigenval_paths = [eigenval_paths; linspace(eigval, eigval, steps)];
        end
    end
    
    for s = 1:steps
        lst = eigvals;
        for idx = 1:size(c_pairs,2)
            pair = c_pairs(:,idx);
            lst(pair(1)) = eigenval_paths(idx,s);
            lst(pair(2)) = 1.0/eigenval_paths(idx,s);
        end
        diag_s = diag(lst);
        path(s,:,:) = (B * diag_s) / B;
    end
    path_out = path;
end

function [c_pairs, eigvals, B] = find_complex_hyperbolic_tuple(S)
    % Returns a tuple of complex eigenvalues of a symplectic matrix S 
    % that do not lie on the unit circle. Eigenvalues that are closer 
    % than eps to the unit circle are discarded.
    % 
    % Parameters
    % ----------
    % S : array
    % Array with shape (2n,2n)
    % 
    % Returns
    % -------
    % Tuple consisting of a list with the indices of complex hyperbolic eigenvalues,
    % a array with the eigenvalues, and 
    % a array with the basis change (from diagonalization)
    % 
    eps_tol = 1e-5;

    [B, eigvals] = eig(S);
    eigvals = diag(eigvals);
    c_pairs = [];
    for i = 1:length(eigvals)
        eigval = eigvals(i);
        if abs(real(eigval) - eigval) < eps_tol || abs(abs(eigval)-1) < eps_tol
            continue
        end
        c_pairs = [c_pairs; i];

%       inverse
        [~, j] = min( abs(eigval-1./eigvals(i+1:end)) );
        inv_idx = i + j;
        c_pairs = [c_pairs; inv_idx];
        % for j = (i+1):length(eigvals)
        %     if abs( 1 / eigvals(j) - eigval ) < eps_tol
        %         c_pairs = [c_pairs; j];
        %         break
        %     end
        % end

%       conjugate
        [~, j] = min( abs( eigval-conj(eigvals(i+1:end)) ) );
        inv_idx = i + j;
        c_pairs = [c_pairs; inv_idx];
        % for j = (i+1):length(eigvals)
        %     if abs( conj(eigvals(j)) - eigval ) < eps_tol
        %         c_pairs = [c_pairs; j];
        %         break
        %     end
        % end
        
%       conjugate inverse
        [~, j] = min( abs( eigval-1./conj(eigvals(i+1:end)) ) );
        inv_idx = i + j;
        c_pairs = [c_pairs; inv_idx];
        % for j = (i+1):length(eigvals)
        %     if abs( 1 / conj(eigvals(j)) - eigval ) < eps_tol
        %         c_pairs = [c_pairs; j];
        %         break
        %     end
        % end

%       pair_found = True
        break
    end 
end


function path_out = path_removing_cpx_hyp_tuple(S, steps)
    % Returns path of symplectic matrices starting at S and ending at S', 
    % where S' has no complex hyperbolic eigenvalues left
    % 
    % Parameters
    % ----------
    % S : array
    % Array with shape (2n,2n)
    % steps : int
    % 
    % Returns
    % -------
    % path : array with shape (steps, 2n, 2n)
    % 
    
    if nargin < 2
      steps = 200;
    end
    
    [tuple, eigvals, B] = find_complex_hyperbolic_tuple(S);
    path = zeros(steps, 4, 4);
    if length(tuple) < 4
        for i = 1:steps
            path(i,:,:) = S;
        end
        path_out = path;
        return
    end
    eigval = eigvals(tuple(1));
    norm_val = norm(eigval);
    angle = atan2(imag(eigval), real(eigval));
    Binv = inv(B);
    angle_path = exp(linspace(angle * 1j, pi * 1j, steps) );
    norm_path = linspace(norm_val, 1.0, steps);
    for s = 1:steps
        lst = eigvals;
        lst(tuple(1)) = norm_path(s) * angle_path(s);
        lst(tuple(2)) = 1 / (norm_path(s) * angle_path(s));
        lst(tuple(3)) = norm_path(s) * conj(angle_path(s) );
        lst(tuple(4)) = 1 / ( norm_path(s) * conj(angle_path(s)) );
        diag_s = diag(lst);
        path(s,:,:) = B * diag_s * Binv;
    end
    path_out = real(path);
end


function R = rotate4(phi)
    c = cos(phi);
    s = sin(phi);
    R = [c, -s, 0, 0;
         s, c, 0, 0;
         0, 0, c, -s;
         0, 0, s, c];
end

% function P = change_order_mat(B, target_order)
%     P = zeros(4);
%     for i=1:4
%         P(i, target_order(i)) = 1.0;
%     end
% end

function [ordering, sorted_vals] = ordering_values(vals)
    sorted_vals = sort(vals, 'descend');
    ordering = [];
    for i=1:length(vals)
        val = sorted_vals(i);
        for j = 1:length(vals)
            if vals(j) == val
                ordering = [ordering j];
            end
        end
    end
end

function path_out = path_removing_positive_hyperbolic_tuple(S, steps)

    if nargin < 2
      steps = 200;
    end
    
    eps_tol = 1e-5;
    
    [c_pairs, eigvals, B] = find_real_hyperbolic_pairs(S);
    % eigvals(abs(imag(eigvals))<eps_tol) = real(eigvals(abs(imag(eigvals))<eps_tol));
    if size(c_pairs,2) < 2 || min(real(eigvals) ) < 0
        path_out = path_removing_cpx_hyp_tuple(S, steps);
        return
    end
    [ordering, eigvals] = ordering_values(eigvals);
    O = zeros(4);
    O(ordering(1),1) = 1;
    O(ordering(4),3) = 1;
    O(ordering(2),2) = 1;
    O(ordering(3),4) = 1;
    BO = B * O;
    BO(:,1) = -1 / sdot4(BO(:,1), BO(:,3)) * BO(:,1);
    BO(:,2) = -1 / sdot4(BO(:,2), BO(:,4)) * BO(:,2);
    BOi = inv(BO);
    ODO = O \ (diag(eigvals) * O);
    ts = linspace(0.0, 1.0, steps);
    lambda_s = linspace(ODO(1,1), 2.0, steps);
    mu_s = linspace(ODO(2,2), 2.0, steps);

    path1 = zeros(steps,4,4);
    path2 = zeros(length(ts),4,4);
    for i = 1:steps
        ls = lambda_s(i); 
        ms = mu_s(i);
        path1(i,:,:) = BO * diag([ls,ms,1/ls,1/ms]) * BOi;
    end
    for i = 1:length(ts)
        s = ts(i);
        path2(i,:,:) = BO * rotate4(0.1*s) * diag([2.0,2.0,0.5,0.5]) * BOi ;
    end
    path3 = path_removing_cpx_hyp_tuple(reshape(path2(end,:,:),4,4), steps);
    path_out = cat(1,path1, path2, path3);
end
    
function path_out = reduce_symplectic_matrix(S, steps)
    if nargin < 2
      steps = 200;
    end
    
    % out_eigs = eig(reshape(S,4,4))
    path2H = path_removing_positive_hyperbolic_tuple(reshape(S,4,4), steps);
    % out_mat = reshape(path2H(end,:,:),4,4)
    % out_eigs = eig(out_mat)
    pathE = path_removing_elliptic_pairs(reshape(path2H(end,:,:),4,4), steps);
    % out_mat = reshape(pathE(end,:,:),4,4)
    % out_eigs = eig(out_mat)
    pathH = path_removing_negative_hyperbolic_pairs(reshape(pathE(end,:,:),4,4), steps);
    % out_mat = reshape(pathH(end,:,:),4,4)
    % out_eigs = eig(out_mat)
    path_out = cat(1,path2H, pathE, pathH);
end



function Ks = interpolate_K(K, dim, steps)
    U = KtoU(K);
    [basis, eigU] = eig(U);
    eigU = diag(eigU);
    basis_inv = inv(basis);
    angles = log(eigU);
    target_angles = zeros(size(angles));
    diags = zeros(length(angles),steps);
    for i = 1:length(angles)
        diags(i,:) = linspace(angles(i), target_angles(i), steps);
    end
    Ks = zeros(steps, dim, dim);
    for i = 1:steps
        Us = diag(exp(diags(:,i)));
        tmp = basis * Us * basis_inv;
        Ks(i,:,:) = UtoK(tmp);
    end
end


function As = interpolate_A(A, dim, steps)
    As = zeros(steps, dim, dim);
    n = dim / 2;
    A_target = eye(n);
    A11_lin = zeros(steps,n,n);
    for i = 1:n
        for j = 1:n
            A11_lin(:,i,j) = linspace(A(i,j), A_target(i,j), steps);
        end
    end
    for i = 1:steps
        As(i,1:n,1:n) = A11_lin(i,:,:);
        As(i,n+1:dim,n+1:dim) = inv(reshape(A11_lin(i,:,:),n,n));
    end
end

function Ns = interpolate_N(N, dim, steps)
    n = dim / 2;
    Ns = zeros(steps, dim, dim);
    N12_lin = zeros(floor(steps/2),n,n);
    for i = 1:n
        for j = 1:n
            N12_lin(:,i,j) = linspace(N(i,j+n), 0.0, floor(steps/2));
        end
    end
    for i = 1:floor(steps/2)
        Ns(i,1:n,1:n) = N(1:n,1:n);
        Ns(i,n+1:dim,n+1:dim) = N(n+1:dim,n+1:dim);
        Ns(i,1:n,n+1:dim) = N12_lin(i,:,:);
    end

    N_lin = zeros(steps-floor(steps/2),dim,dim);
    target_N = eye(dim);
    for i = 1:dim
        for j = 1:dim
            N_lin(:,i,j) = linspace(Ns(floor(steps/2),i,j), target_N(i,j), steps-floor(steps/2));
        end
    end
    for i  = (floor(steps/2)+1):steps
        Ns(i,:,:) = N_lin(i - floor(steps/2),:,:);
    end
end



function path_out = path_hyperbolic_basepoint(S, steps)
    if nargin < 2
      steps = 200;
    end
    Omega4 = get_omega0(4);

    [c_pairs, eigvals, B] = find_real_hyperbolic_pairs(S);
    if size(c_pairs,2) ~= 1
        path_out = zeros(1, 4, 4);
        path_out(1,:,:) = S;
        return
    end
    pair = c_pairs(:,1);
    B1 = B;
    B1(:,pair(1)) = B1(:,pair(1)) / (-dot(B(:,pair(1)), Omega4 * B(:,pair(2)) ));
    other_pair = [1, 2, 3, 4];
    other_pair = other_pair(~ismember(other_pair, pair));
    
    B1(:,other_pair(1)) = B1(:,other_pair(1)) / (-dot(B(:,other_pair(1)), Omega4 * B(:,other_pair(2)) ));
    
    B1 = real(B1);
    first_pair = other_pair;
    second_pair = pair;
    if ismember(1, pair)
        first_pair = pair;
        second_pair = other_pair;
    end
    
    idx_order = [1 second_pair(1) first_pair(2) second_pair(2)];
    B2 = zeros(4);
    R = zeros(4);
    for i = 1:4
        B2(:,i) = B1(:, idx_order(i));
        R(i,idx_order(i)) = 1.0;
    end
    Rinv = inv(R);
    B1 = B1 * R ;
    [K, A, N] = iwasawa(B1);
    Kpath = interpolate_K(K, length(K), steps);
    Apath = interpolate_A(A, length(A), steps);
    Npath = interpolate_N(N, length(N), steps);

    Bpath = zeros(steps, 4, 4);
    for i = 1:steps
        Bpath(i,:,:) = reshape(Kpath(i,:,:),4,4) * reshape(Apath(i,:,:),4,4) * reshape(Npath(i,:,:),4,4);
    end
    Spath = zeros(steps, 4, 4);
    for i = 1:steps
        Bi = reshape(Bpath(i,:,:),4,4) * Rinv;
        % Bi_inv = inv(Bi);
        Spath(i,:,:) = real( (Bi * diag(eigvals)) / Bi );
    end
    nf = reshape(Spath(end,:,:),4,4); % normal form
    [max_eigval, argmax]=max(reshape(nf',16,1)); %To row major, like Python
    idx0 = mod(argmax-1, 4)+1;
    if idx0 + 2 < 5
        idx1 = idx0 + 2;
    else
        idx1 = idx0 - 2;
    end
    rot_steps = floor(steps / 10);
    eigenval_path = linspace(max_eigval, 2.0, rot_steps);
    Rpath = zeros(rot_steps, 4, 4);
    for i = 1:rot_steps
        Rpath(i,:,:) = nf;
        Rpath(i,idx0, idx0) = eigenval_path(i);
        Rpath(i,idx1, idx1) = 1.0 / eigenval_path(i);
    end
    path_out = real( cat(1,Spath, Rpath) );
end

function [ang_idx, angle_jumps] = compute_angular_index(cdets, T)
    angle_jumps = imag(log( cdets(2:T)./cdets(1:T-1) ));
    fprintf("\nSum of angle jumps divided by pi (should be very near integer): %.10f \n", sum(angle_jumps) /pi)
    ang_idx = round(sum(angle_jumps) /pi);
end

function out_vec = check_symplecticity(s_path)
    Omega4 = get_omega0(4);
    out_vec = zeros(1,size(s_path,1));
    for i = 1:size(s_path,1)
        A = reshape(s_path(i,:,:),4,4);
        out_vec(i) = sum(sum(abs(A' * Omega4 * A - Omega4)));
    end
end

function out_vec = check_maslov_intersections(s_path)
    out_vec = zeros(1,size(s_path,1) - 1);
    for i = 1:size(s_path,1) - 1
        out_vec(i) = maslov_func(reshape(s_path(i+1,:,:),4,4)) * maslov_func(reshape(s_path(i,:,:),4,4));
    end
end

function out_vec = check_continuity(s_path)
    out_vec = zeros(1,size(s_path,1)-1);
    for i = 1:size(s_path,1)-1
        out_vec(i) = sum(sum( abs(s_path(i+1,:,:)-s_path(i,:,:)) ));
    end
end

function [CZ, extend_s_path, path_reduced, second_path] = get_index_sympl_path(s_path, steps, error_report)
    
    if nargin < 2
      steps = 20000;
      error_report = false;
    elseif nargin < 3
      error_report = false;
    end

    path_reduced = reduce_symplectic_matrix(reshape(s_path(end,:,:),4,4), steps);
    S_reduced = real(path_reduced(end,:,:));
    second_path = path_hyperbolic_basepoint(reshape(S_reduced,4,4), steps);
    concat_path = cat(1,path_reduced, second_path);
    extend_s_path = cat(1, s_path, concat_path );
    cdets = retract_sympl_path(extend_s_path);
    [CZ, angle_jumps] = compute_angular_index(cdets, length(cdets));
    % For reliability checking
    if error_report
        symplecticity = max( check_symplecticity(extend_s_path));
        extension_sign = min(check_maslov_intersections(concat_path));
        continuity = max(check_continuity(extend_s_path));
        [max_angle_jump, idx_jump] = max(angle_jumps);
        fprintf("Max angle jump at %d, %d; jump size = %.12f \n", idx_jump, length(angle_jumps), max_angle_jump);
        if max_angle_jump > pi
            fprintf("Angle jump is too large: increase the number of integration steps\n");
        end
        fprintf("Error in symplecticity was at most %.12f \n", symplecticity)
        if extension_sign > 0
            fprintf("Extention stayed away from Maslov cycle; closest at value %.12f \n", extension_sign);
        else
            fprintf("Extention crossed Maslov cycle; refine number of steps, eps %.12f \n", extension_sign);
        end
        fprintf("Error in continuity was at most %.12f; increase number of steps if this is too large.\n", continuity);
    end
    
end
    

