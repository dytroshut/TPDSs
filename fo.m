function Xmat = fo(X,n)
%n is the dimension of the third mode
[rn,c]=size(X);
r = rn/n;
Xmat = zeros(r,c,n);
for i = 1:n
    Xmat(:,:,i) = X((i-1)*r+1:i*r,:);
end