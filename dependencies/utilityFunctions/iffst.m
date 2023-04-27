function y = iffst(alpha,imSize,Sk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFST Vectorized inverse fast finite shearlet transform
% IN    alpha:  vectorized shearlet coefficients
%       imSize: 2-element array with image dimensions (T,M)
%       Sk:     dictionary bases
% -----------------------------------------------------------------------------------
% OUT   y:      vectorized 2D image
% -----------------------------------------------------------------------------------
% AUTHOR: E. Zea (zea@kth.se)
% DATE: 2019-03-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha    = reshape(alpha,imSize(1),imSize(2),size(Sk,3));   % reshape to 3D array...
sum_alpha = sum(fftshift(fftshift(fft2(alpha),1),2).*Sk,3); % sum up decompositions
y = real(ifft2(ifftshift(sum_alpha))); % recover RIR image
y = y(:); % vectorize output...
end