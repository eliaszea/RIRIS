function [image_ext,imSize,M0,LAMBDA,origAprtr] = ...
                        extRIRimage(image_under,lambda,ranPos)
% EXTRIRIMAGE Spatial extrapolation of RIR image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IN    image_under: under-sampled RIR image
%       lambda:      Hadamard masking operator
%       ranPos:      random mic array positions
% -----------------------------------------------------------------------------------
% OUT   image_ext:   extrapolated RIR image
%       imSize:      2-element array with 
%       M0:          original no. mic positions
%       LAMBDA:      updated masking matrix (includes extrapolation)
%       origAprtr:   original aperture (indices)
% -----------------------------------------------------------------------------------
% Note: by default M is extrapolated to the next power of 2 w.r.t M0
% -----------------------------------------------------------------------------------
% AUTHOR: E. Zea (zea@kth.se)
% DATE: 2019-03-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
padFlag = 1; % lpbp extrap. (default true)
T = size(image_under,1); % fetch no. time samples
if padFlag
    M0 = size(image_under,2); % save original no. mic positions
    M = 2^nextpow2(M0); % update no. mic positions (next power of 2->fft-efficient!)
    origAprtr = (M-M0)/2+1:M-(M-M0)/2; % save original array aperture
    image_ext = zeros(T,M); % preallocate extrapolated image
    for tt = 1:T  % execute linear predictive border padding
        image_ext(tt,:) = lpbp1D(image_under(tt,:),M,2);
    end
    lambda_ext = ones(T*(M-M0)/2,1); % extrapolated data is not masked!
    lambda = [lambda_ext; ...
                lambda(:); ...
                    lambda_ext]; % concatenate before and after lambda
else
    origAprtr = 1:M; % array aperture does not change
    image_ext = image_under; % if no padding just assign same input image
end
imSize = [T,M]; % update image size
LAMBDA = sparse(1:numel(lambda(:)),...
                1:numel(lambda(:)),lambda(:)); % compute sparse masking mtx
zerosUnder = setdiff(1:M,ranPos+(M-M0)/2);   % set of under-sampled pixels
set = []; 
for ii = 1:numel(zerosUnder) % update set of rows to be removed
    set = [set; (zerosUnder(ii)-1)*T+1:zerosUnder(ii)*T];
end
LAMBDA(set,:) = [];
end