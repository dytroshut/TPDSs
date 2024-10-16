%% Check if A = X1 * X^(+)
n = 2;
h = 2;
m = 7;
r = 2^m;
l = 10;
%X0 = {x(0), x(1),x(2),...,x(l)}
x0 = randn(n,h,r);
%initialize x0
X = zeros(n,l*h,r);
X1 = X;
%each x(i) is n*h*r and we have l x(i)'s
X(:,1:h,:) = x0;
A = randn(n,n,r)/50;
bA = bcir(A);
%initialize A
for i = 1:l-1
    v = X(:,(i-1)*h+1:i*h,:);
    X(:,i*h+1:(i+1)*h,:) = fo(bA*unf(v),r);
    X1(:,(i-1)*h+1:i*h,:) = X(:,i*h+1:(i+1)*h,:);
    %to calculate x(i), we use x(i) = A * x(i-1) = 
    % fold(bcirc(A) unfold(x(i-1))
end
v = X(:,(l-1)*h+1:l*h,:);
X1(:,(l-1)*h+1:l*h,:) = fo(bA*unf(v),r);
Xp = unbir(pinv(bcir(X)),r);
B = tprod(X,Xp);
C = tprod(X1,Xp);
norm(bcir(A) - bcir(tprod(X1,Xp)))
% A should be the same as tprod(X1,Xp);
%% Part B stability 
le = 13;
n = 2;
h = 2;
l = 10;
itel = 10;
t1 = zeros(itel,le);
t2 = t1;
dif = zeros(itel);
for m = 1:le
for ite = 1:itel
r = 2^m;
sig = zeros(n*r,itel);
sigma = zeros(n*r,itel);
%X0 = {x(0), x(1),x(2),...,x(l)}
x0 = randn(n,h,r);
%initialize x0
X = zeros(n,l*h,r);
X1 = X;
%each x(i) is n*h*r and we have l x(i)'s
X(:,1:h,:) = x0;
A = randn(n,n,r)/50;
bA = bcir(A);
%initialize A
for i = 1:l-1
    v = X(:,(i-1)*h+1:i*h,:);
    X(:,i*h+1:(i+1)*h,:) = fo(bA*unf(v),r);
    X1(:,(i-1)*h+1:i*h,:) = X(:,i*h+1:(i+1)*h,:);
    %to calculate x(i), we use x(i) = A * x(i-1) = 
    % fold(bcirc(A) unfold(x(i-1))
end
v = X(:,(l-1)*h+1:l*h,:);
X1(:,(l-1)*h+1:l*h,:) = fo(bA*unf(v),r);
tic
sigma(:,ite) = sort(eig(bA));
%calculate the eigenvalues directly from the circulant matrix
t1(ite,m) = toc
tic
fA = fft(A,[],3);
%calculate the sum of the rank of block matrices from fft(bcirc(X0))
for i = 1:r
    sig(2*i-1:2*i) = eig(fA(:,:,i));
    %calculate eigenvalues of each small diagonal matrices
end
sig = sort(sig);
t2(ite,m) = toc
end
end
T1 = mean(t1);
T2 = mean(t2);
y = 1:le;
R = 2.^y;
loglog(R, T1,R,T2,'r-*')
% do the loglog plot
legend('prop4','cor6')
xlabel('r');
ylabel('time');
title('stability')
tt1 = T1./R.^3;
tt2 = T2./R;
figure
plot(R(4:end), tt1(4:end) ,'r-*')
legend('prop4 time/r^3')
xlabel('r');
ylabel('ratio');

figure
plot(R(4:end), tt2(4:end),'r-*')
legend('cor6 time/r')
xlabel('r');
ylabel('ratio');
T1./R.^3
T2./R
%to show T1 is O(r^3) and T2 is O(r)
%% check if the eigenvalues got from two methods are close to each other 
ite = 3;
n = 2;
r = 2^13;
dif = zeros(1,ite);
for j = 1:ite
sig = zeros(n*r,ite);
sigma = sig;
tic
A = randn(n,n,r)/50;
bA = bcir(A);
sigma(:,j) = sort(eig(bA));

fA = fft(A,[],3);
%calculate the sum of the rank of block matrices from fft(bcirc(X0))
for i = 1:r
    sig(2*i-1:2*i,j) = eig(fA(:,:,i));   
end
sig(:,j) = sort(sig(:,j));
dif(j) = norm(sigma(:,j)-sig(:,j))/norm(sigma(:,j));
toc
end
dif