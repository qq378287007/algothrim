function  beta = myLMfit(X, Y, model, beta, maxiter)
% Levenberg-Marquardt algorithm for nonlinear regression

% Set up convergence tolerances from options.
betatol = 1.00E-08;
rtol = 1.00E-08;

% Set the iteration step
sqrteps = sqrt(eps);

% Set initial weight for LM algorithm.
lambda = 0.01;

% treatment for nans
yfit = model(beta, X);%求值
r = Y(:) - yfit(:);%求差
sse = r' * r;%误差平方和

for iter = 1 : maxiter
    betaold = beta;
    sseold = sse;
    
    % Compute a finite difference approximation to the Jacobian
    J = getjacobian(X, model, beta);
    
    % Levenberg-Marquardt step: inv(J'*J+lambda*D)*J'*r
    diagJtJ = sum(J .^ 2, 1);
    Jplus = [J; diag(sqrt(lambda * diagJtJ))];
    rplus = [r; zeros(length(beta), 1)];
    step = Jplus \ rplus;% Jplus * step = rplus，求step
    %step = pinv(Jplus) * rplus;
    %step = pinv(J' * J + lambda * diagJtJ) * J' * r;
    %step = (Jplus' * Jplus) ^ -1 * Jplus' * rplus;
    beta(:) = beta(:) + step;
    
    % Evaluate the fitted values at the new coefficients and
    % compute the residuals and the SSE.
    yfit = model(beta, X);
    r = Y(:) - yfit(:);
    sse = r' * r;
    
    % If the LM step decreased the SSE, decrease lambda to downweight the
    % steepest descent direction.  Prevent underflowing to zero after many
    % successful steps; smaller than eps is effectively zero anyway.
    if sse < sseold
        lambda = max(0.1*lambda, eps);
        
        % If the LM step increased the SSE, repeatedly increase lambda to
        % upweight the steepest descent direction and decrease the step size
        % until we get a step that does decrease SSE.
    else
        while sse > sseold
            lambda = 10 * lambda;
            if lambda > 1e16
                return;
            end
            
            %diagJtJ = sum(J .^ 2, 1);
            Jplus = [J; diag(sqrt(lambda * diagJtJ))];
            
            % 最小二乘矩阵解，A * x = b, x = (A' * A) ^ -1 * A' * b
            % Jplus * step = rplus，求step
            step = Jplus \ rplus;
            %step = pinv(J' * J + lambda * diagJtJ) * J' * r;
            %step = pinv(Jplus) * rplus;
            %广义逆矩阵pinv(Jplus)不唯一，(Jplus' * Jplus) ^ -1 * Jplus'
            %step = (Jplus' * Jplus) ^ -1 * Jplus' * rplus;
            beta(:) = betaold(:) + step;
            
            yfit = model(beta, X);
            r = Y(:) - yfit(:);
            sse = r' * r;
        end
    end
    
    % Check step size and change in SSE for convergence.
    %norm, 平方和开方
    if norm(step) < betatol * (sqrteps + norm(beta)) ...
            || abs(sse - sseold) <= rtol * sse
        break;
    end
end
