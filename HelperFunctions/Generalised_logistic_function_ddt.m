function dpdt = Generalised_logistic_function_ddt(t,t0,tau,A)
% https://en.wikipedia.org/wiki/Generalised_logistic_function
    %p = (1+exp(-(t-t0)/tau)).^(-A);
    dpdt = (A*exp(-(t - t0)/tau))./(tau*(exp(-(t - t0)/tau) + 1).^(A + 1));
end



