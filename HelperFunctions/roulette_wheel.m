function index = roulette_wheel(w,N)
%roulette_wheel Summary of this function goes here
%   https://stackoverflow.com/questions/2977497/weighted-random-numbers-in-matlab
    a = 1:numel(w);
    index = a( sum( bsxfun(@ge, rand(N,1), cumsum(w./sum(w))), 2) + 1 );
end

