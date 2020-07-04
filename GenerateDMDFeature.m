function [sorted_modes, sorted_lam, sorted_power, freq] = GenerateDMDFeature(X,Y,dt,nb_mode)
% This function calculates DMD and returns DMD mode, lambda, frequency and power.
% Inputs:
% X -> X_1: [x_1, x_2, ..., x_{m-1}] data matrix
% Y -> X_2: [x_2, x_3, ..., x_m] data matrix
% dt -> 1/fs Time per sample
% nb_mode -> The number of DMD modes you want to calculate
% Outputs:
% sorted_modes -> DMD modes sorted by frequency
% sorted_lam -> eigenvalues sorted by frequency
% sorted_power -> powers sorted by frequency
% freq -> frequency

%[U,S,V] = svd(X, 'econ');
[U,S,V] = svds(X, nb_mode);
M = Y*V/S;
A = U'*M;

% Ahat is normalized A.
Ahat = S^(-1/2) * A * S^(1/2);
[What, Dhat,zhat] = eig(Ahat);
W = S^(1/2) * What;

v = (M*W)/Dhat; % exact DMD

lam = diag(Dhat);

omega = log(lam) / dt;
freq = abs(imag(omega) / (2*pi));
Phi = v * diag(lam);
Power = diag(Phi'*Phi);

% sort DMD modes
temp = horzcat(freq,v');
temp = sortrows(temp,1);
temp(:,1) = [];
sorted_modes = temp;
sorted_modes = sorted_modes';
clear temp

% sort lambdas
temp = horzcat(freq, lam);
temp = sortrows(temp,1);
sorted_lam = temp(:,2);
clear temp

% sort powers
temp = horzcat(freq,Power);
temp = sortrows(temp,1);
sorted_power = temp(:,2);

% sort frequencies
freq = sortrows(freq);