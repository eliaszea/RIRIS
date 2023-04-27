function alpha = ffst(y,imSize,Sk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FFST Vectorized fast finite shearlet transform
% IN    y:      vectorized 2D image
%       imSize: 2-element array with image dimensions (T,M)
%       Sk:     dictionary bases
% -----------------------------------------------------------------------------------
% OUT   alpha:  vectorized shearlet coefficients
% -----------------------------------------------------------------------------------
% AUTHOR: E. Zea (zea@kth.se)
% DATE: 2019-03-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = reshape(y,imSize(1),imSize(2)); % reshape to image input...
uST = Sk.*repmat(fftshift(fft2(y)),[1,1,size(Sk,3)]); % 2D Fourier multiplications
alpha = real(ifft2(ifftshift(ifftshift(uST,1),2))); % space-time decomposition
alpha = alpha(:); % vectorize output...
end