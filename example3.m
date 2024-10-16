%tensor last example
n = 2;
h = 2;
l = 10;
ite =10;
%iteration 10 times 
le = 8;
%the third mode dimension goes from 2 to 2^le
ras = zeros(ite,le); ras2 = zeros(ite,le);
%ras and ras2 saves the total rank and we can compare
%them with nr to see if the system is informative
t1 = zeros(ite,le); %t2 = ra;
%t3 = t1;
t2 = t1;
%t1 and t2 store the time spent to calculate the total rank
%t1(i,j) means ith experiment for r = 2^j (third mode dimension)
for m = 1:le
for q = 1:ite
r = 2^m;
%X0 = {x(0), x(1),x(2),...,x(l)}
x0 = randn(n,h,r);
%initialize x0
X = zeros(n,l*h,r);
%each x(i) is n*h*r and we have l x(i)'s
X1 = X;
X(:,1:h,:) = x0;
A = randn(n,n,r)/50;
%initialize A
bA = bcir(A);
for i = 1:l-1
    v = X(:,(i-1)*h+1:i*h,:);
    X(:,i*h+1:(i+1)*h,:) = fo(bA*unf(v),r);
    X1(:,(i-1)*h+1:i*h,:) = X(:,i*h+1:(i+1)*h,:);
    %to calculate x(i), we use x(i) = A * x(i-1) = 
    % fold(bcirc(A) unfold(x(i-1))
end
v = X(:,(l-1)*h+1:l*h,:);
X1(:,(l-1)*h+1:l*h,:) = fo(bA*unf(v),r);

% Define 'a' as a symbolic variable
syms a;
% Convert X1 and X to symbolic
X1_sym = sym(X1);
X_sym = sym(X);
% Create the matrix expression X1 - a * X
matrix_expr = X1_sym - a * X_sym;
% Calculate the symbolic rank of the matrix expression
%rank_expr = rank(bcir(matrix_expr));
% Substitute a specific value for 'a', for example a = 1
%a_value = 1;
%rank_at_a_1 = subs(rank_expr, a, a_value);
tic
ras(q,m) = rank(bcir(matrix_expr));
%calculate rank of bcirc(X0)
t1(q,m) = toc


%tic
%Xfc = ffcir(X);

%for i = 1:r
%    Ti = Xfc((i-1)*n+1:i*n,(i-1)*l*h+1:i*l*h);
%    ras(m) = ras(m) + rank(Ti);
%end
%t2(m) = toc
tic
fX1 = fft(X1,[],3);
fX = fft(X,[],3);
fX1_sym = sym(fX1);
fX_sym = sym(fX);
fmatrix_expr = fX1_sym - a * fX_sym;
%calculate the sum of the rank of block matrices from fft(bcirc(X0))
for i = 1:r
    Ti2 = fmatrix_expr(:,:,i);
    ras2(q,m) = ras2(q,m) + rank(Ti2);
end
t2(q,m) = toc

%tic
%fX = symbolic_fft(matrix_expr,3);
%fX = symbolic_fft_nd(matrix_expr,3);
%calculate the sum of the rank of block matrices from fft(bcirc(X0))
%for i = 1:r
%    Ti2 = fX(:,:,i);
%    ras2(q,m) = ras2(q,m) + rank(Ti2);
%end
%t3(q,m) = toc
%T3(m) = T3(m)+t3(q,m);
y = 1:le;
R = 2.^y;
end
end
T1 = mean(t1);
T2 = mean(t2);
y = 1:le;
R = 2.^y;
loglog(R, T1,R,T2,'r-*')
% do the loglog plot
legend('prop7','cor9')
xlabel('r');
ylabel('time');
title('stability')
tt1 = T1./R.^3;
tt2 = T2./R;
figure
plot(R, tt1 ,'r-*')
legend('prop7 time/r^3')
xlabel('r');
ylabel('ratio');

figure
plot(R, tt2,'r-*')
legend('cor9 time/r')
xlabel('r');
ylabel('ratio');
T1./R.^3
T2./R
%to show T1 is O(r^3) and T2 is O(r)
