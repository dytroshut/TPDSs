function Xmat = unf(X)
[r,c,h]=size(X);
Xmat = zeros(r*h,c);
for i = 1:h
    Xmat((i-1)*r+1:i*r,:) = X(:,:,i);
end