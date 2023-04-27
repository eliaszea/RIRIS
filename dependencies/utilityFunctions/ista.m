function alpha = ista(y_hat,LAMBDA,E,nu_max,imSize,Sk,beta,epsilon)
% ISTA Iterative soft thresholding algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IN    y_hat:   observation vector
%       LAMBDA:  masking matrix
%       E:       expansion matrix
%       nu_max:  max no. iterations
%       imSize:  2-element array with image dimensions (T,M)
%       Sk:      dictionary bases
%       beta:    threshold choice
%       epsilon: estimated noise level
% -----------------------------------------------------------------------------------
% OUT   alpha:   thresholded shearlet coefficients
% -----------------------------------------------------------------------------------
% AUTHOR: E. Zea (zea@kth.se)
% DATE: 2019-03-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nStarting thresholding iterations...\n');
alpha = E'*ffst(LAMBDA'*y_hat,imSize,Sk); % initialize solution vector
its = 0; % init iteration counter
C = linspace(1,epsilon/beta,nu_max+1);  % decreasing threshold function
zeta0 = beta*norm(alpha,Inf);           % first zeta
res_norm = norm(LAMBDA*iffst(E*alpha,imSize,Sk)-y_hat); % initial residual norm
while ( res_norm > epsilon ) && ( its < nu_max )
    % update iteration count
    its = its + 1; 
    % update solution 
    alpha = wthresh(alpha + E'*ffst( LAMBDA'*(y_hat - ...
                    LAMBDA*iffst(E*alpha,imSize,Sk)),imSize,Sk ),...
                    's',C(its)*zeta0); 
	% update residual norm
    res_norm  = norm(LAMBDA*iffst(E*alpha,imSize,Sk)-y_hat); 
    disp(['Iteration no. ' num2str(its) '/' num2str(nu_max)]);
end
fprintf('Interpolating image from thresholded coefficients...\n');
fprintf('\n---------------- INTERPOLATION DONE! ----------------\n');
end