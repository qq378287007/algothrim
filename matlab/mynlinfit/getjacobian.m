function J = getjacobian(X, model, beta)
    function yplus = call_model_nested(betaNew)
        yplus = model(betaNew, X);
    end
J = statjacobian(@call_model_nested, beta);
end
