function [lambda,ranPos,N] = selectionMtrx(u,imSize)
% SELECTIONMTRX Generate masking (selection) matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IN    u:          under-sampling factor
%       imSize:     2-element array with image dimensions (T, M)
% -----------------------------------------------------------------------------------
% OUT   lambda:     Hadamard (point-wise) mask 
%       ranPos:     random mic positions
%       N:          no. mic positions (post-under-sampling)
% -----------------------------------------------------------------------------------
% AUTHOR: E. Zea (zea@kth.se)
% DATE: 2019-03-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jitter = u; % jittering (default = u)
rng(35111); % set random number generator (changes array topology!)
if mod(u,2)
    f0 = (1-u)/2;
else
    f0 = -u/2;
end
mu = round(imSize(2)/u); % apply under-sampling factor to spatial dimension
f = f0 + u*[1:mu]';      % initialize deterministic and stochastic parts
epsilon = randi([-floor((jitter-1)/2) floor((jitter-1)/2)],mu,1); 
ranPos = f + epsilon;    % random microphone positions to be removed
N = numel(ranPos);       % no. of random measurement samples
% masking operation in Hadamard form
lambda = zeros(imSize(1),imSize(2));
for rr = 1:numel(ranPos)        
    lambda(:,ranPos(rr)) = ones(imSize(1),1);
end
end