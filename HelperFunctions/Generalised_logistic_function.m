function p = Generalised_logistic_function(t,t0,tau,A)
% https://en.wikipedia.org/wiki/Generalised_logistic_function
    p = (1+exp(-(t-t0)/tau)).^(-A);
end



