function signalOut = lpbp1D( signalIn, Ntot, order )
%LPBP1D 1D Linear Predictive Border Padding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function extrapolates linear microphone array data by means of
%   designing and applying filters with AR prediction coefficients. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IN    signalIn:   1D signal in a vector
%       Ntot:       no. of total samples (after extrapolation)
%       order:      order of the filter ("no. of Fourier peaks")
% -----------------------------------------------------------------------------------
%   OBS: The difference in samples between Ntot and the length of 
%   signalIn MUST be an even number!
%
% OUT   signalOut:  extrapolated signal
% -----------------------------------------------------------------------------------
% EXAMPLE OF USAGE
%   Apply a 4-th order LPBP extrapolation to the next power-of-two samples
%       
%       >> p_ext = lpbp1D( p, 2^(nextpow(length(p)), 4 );
% -----------------------------------------------------------------------------------
% REFERENCES 
%   [1] R. Scholte et al. Truncated aperture extrapolation for Fourier-based 
%       near-field acoustic holography by means of border-padding, JASA 125(6) 
%       pp. 3844-3854 (2009)
% -----------------------------------------------------------------------------------
% AUTHOR: E. Zea (zea@kth.se)
% DATE: 2019-03-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% difference in samples between in and out
dN = 0.5*(Ntot-length(signalIn)); 
if mod(dN,2) ~= 0
    error('The extrapolated number of samples must be even!');
end
signalOut = zeros( Ntot, 1 );           % preallocate output signal
signalOut(dN+1:Ntot-dN,:) = signalIn;   % allocate input signal in the middle
signal_extrapolated_flip = flipud(signalOut); % needs flip first!
% backwards extrapolation
a_back = arburg(signal_extrapolated_flip(dN+1:Ntot-dN,1),order);
Z_back = filtic(1,a_back,signal_extrapolated_flip(Ntot-dN-(0:(order-1)),1));
signalOut(1:dN,1) = fliplr(filter(1,a_back,zeros(1,dN),Z_back));
% forward extrapolation
a_forw = arburg(signalOut(dN+1:Ntot-dN,1),order);
Z_forw = filtic(1,a_forw,signalOut(Ntot-dN-(0:(order-1)),1));
signalOut(Ntot-dN+1:Ntot,1) = filter(1,a_forw,zeros(1,dN),Z_forw);
end