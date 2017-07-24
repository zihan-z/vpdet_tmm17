function [Matrix,Index] = npermutek(N,K)
%NPERMUTEK(N,K) Permutations of elements of N, K at a time, with repetition.
% MAT = NPERMUTEK(N,K) returns all possible samplings from vector N of  
% type: ordered sample with replacement.  
% MAT has size (length(N)^K)-by-K, where K must be a scalar.
% [MAT, IDX] = NPERMUTEK(N,K) also returns IDX such that MAT = N(IDX).
%
% Example:
%         MAT = npermutek([2 4 5],2)
%
%  MAT =
% 
%       2     2
%       2     4
%       2     5
%       4     2
%       4     4
%       4     5
%       5     2
%       5     4
%       5     5
%
%
% See also perms, nchoosek
%
% Also on the web:
% http://mathworld.wolfram.com/BallPicking.html
% See the section on Enumerative combinatorics below: 
% http://en.wikipedia.org/wiki/Permutations_and_combinations
% Author:  Matt Fig
% Contact:  popkenai@yahoo.com

if nargin ~= 2
    error('NPERMUTEK requires two arguments. See help.')
end

if isempty(N) || K == 0,
   Matrix = [];  
   Index = Matrix;
   return
elseif floor(K) ~= K || K<0 || ~isreal(K) || numel(K)~=1 
    error('Second argument should be a real positive integer. See help.')
end

if K==1
    Matrix = N(:); % This one is easy to calculate.
    Index = (1:length(N))';
    return
end

lgth = numel(N); % Used in calculating the Index.
Index = zeros(lgth^K,K); % Preallocation.

Index(:,K) = 1; % We don't need to do these in loop.
Index(1:lgth^(K-1):end,1) = 1;
Index(lgth+1:lgth:end,K) = 1-lgth;

for ii = K-1:-1:2
    Index(1:lgth^(ii-1):end,K-ii+1) = 1; % Poke in ones.
    Index(lgth^ii+1:lgth^ii:end,K-ii+1) = 1-lgth; % Poke in drops.
end

Index = cumsum(Index,1); % Create the Index matrix.
Matrix = N(Index);  % This takes most of the runtime.



