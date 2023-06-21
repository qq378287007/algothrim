function beta = mynlinfit(X, Y, model, beta)
maxiter = 200;
beta = myLMfit(X, Y, model, beta, maxiter);
