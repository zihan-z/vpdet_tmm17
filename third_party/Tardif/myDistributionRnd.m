%desctrete distribution with probabilities of k - p[1, k]
function res = myDistributionRnd( p )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

k = size (p, 2);
l = zeros(1, 1:k); r = zeros (1, 1:k);
l(1, 1) = 0;
r(1, 1) = p(1, 1);
for i = 2:k
    l(1, i) = p(1, i-1);
    r(1, i) = l(1, i) + p(1, i);
end

rndn = rand(1);

for i = 1:k
    if (rndn < r(1, i) && rndn >= l(1, i)
        res = i;
    end
end    

end

