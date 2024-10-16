function Xcir = bcir(X)
[r,c,n] = size(X);
Xcir = zeros(r*n,c*n);
v = 1:n;
V = zeros(n);
for i = 0:n-1
    V(i+1, :) = circshift(v, [0, i]); % Shift vector to create rows
end
V = V';
for i = 1:n
    for j = 1:n
        h = V(i,j);
        Xcir((i-1)*r+1:i*r,(j-1)*c+1:j*c) = X(:,:,h);
    end
end