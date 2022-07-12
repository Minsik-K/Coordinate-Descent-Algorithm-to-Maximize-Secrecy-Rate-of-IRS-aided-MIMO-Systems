function [Covariance_matrix] = MIMO_Capacity_geneig(w, H, R, A, lambda)
% 
% [1] D. Park, "Weighted sum rate maximization of MIMO broadcast and
% interference channels with confidential messages," IEEE Trans. Wireless
% Commun., vol. 15, no. 3, pp. 1742-1753, Mar. 2016.
%
% [2] S.-J. Kim and G. B. Giannakis, "Optimal resource allocation for MIMO
% ad hoc cognitive radio networks," IEEE Trans. Inf. Theory, May 2011.
%
% Lemma 1 in [1] is equivalent to Proposition 1 in [2].
%

% Lemma 1 in [1]
M=size(H,2);
A_inv = inv(sqrtm(A+lambda*eye(M)));
[~, S, V] = svd( inv(sqrtm(R))*H*A_inv, 'econ');
s = real(diag(S));
P = w-1./s.^2;
P = diag(max(P,0));
Covariance_matrix=A_inv*V*P*V'*A_inv;