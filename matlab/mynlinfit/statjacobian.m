function J = statjacobian(func, beta)
%STATJACOBIAN Estimate the Jacobian of a function

% Handle optional arguments, starting with y0 since it will be needed to
% determine the appropriate size for a default groups.

y0 = func(beta);

% When there is only one group, ensure that beta is a row vector so
% that vectoriation works properly. Also ensure that the underlying
% function is called with an input with the original size of beta.
thetaOriginalSize = size(beta);
beta = reshape(beta, 1, []);

func = @(beta) func(reshape(beta, thetaOriginalSize));

% Best practice for forward/backward differences:
DerivStep = eps^(1/3);
nTheta = sqrt(sum(beta(:) .^ 2));
DerivStep2 = DerivStep * (nTheta + (nTheta == 0.0));

[~, numParams] = size(beta);
J = zeros(numel(y0), numParams);
for i = 1 : numParams
    % Calculate delta(:, i), but remember to set it back to 0 at the end of the loop.
    delta = DerivStep * beta(i);
    if delta == 0.0
        delta = DerivStep2;
    end
    thetaNew = beta;
    thetaNew(i) = thetaNew(i) + delta;
    yplus = func(thetaNew);
    dy = yplus(:) - y0(:);
    J(:, i) = dy ./ delta;
end

