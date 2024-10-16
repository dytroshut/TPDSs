function X = unbir(bX,s)
%given bcirc(X) with dimension ns by ms, and given s
%X is n by m by s
[a,b] = size(bX);
n = a/s; m = b/s;
X = zeros(n,m,s);
for i = 1:s
    X(:,:,i) = bX((i-1)*n+1:i*n,1:m);
end