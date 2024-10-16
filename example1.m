ite = 10;
%iteration 10 times 
le = 12;
%the third mode dimension goes from 2 to 2^le
ra = zeros(1,le);
ras = ra; ras2 = ra;
%ras and ras2 saves the total rank and we can compare
%them with nr to see if the system is informative
t1 = zeros(ite,le); %t2 = ra;
t3 = t1;
%t1 and t3 store the time spent to calculate the total rank
%t1(i,j) means ith experiment for r = 2^j (third mode dimension)
T1 = ra; T3 = ra;
%T1 and T3 are the average time spent
for m = 1:le
for q = 1:ite
n = 2;
h = 2;
r = 2^m;
l = 10;
%X0 = {x(0), x(1),x(2),...,x(l)}
x0 = randn(n,h,r);
%initialize x0
X = zeros(n,l*h,r);
%each x(i) is n*h*r and we have l x(i)'s
X(:,1:h,:) = x0;
A = randn(n,n,r);
%initialize A
for i = 1:l-1
    v = X(:,(i-1)*h+1:i*h,:);
    X(:,i*h+1:(i+1)*h,:) = fo(bcir(A)*unf(v),r);
    %to calculate x(i), we use x(i) = A * x(i-1) = 
    % fold(bcirc(A) unfold(x(i-1))
end
tic
ras(m) = rank(bcir(X));
%calculate rank of bcirc(X0)
t1(q,m) = toc
T1(m) = T1(m)+t1(q,m);


%tic
%Xfc = ffcir(X);

%for i = 1:r
%    Ti = Xfc((i-1)*n+1:i*n,(i-1)*l*h+1:i*l*h);
%    ras(m) = ras(m) + rank(Ti);
%end
%t2(m) = toc

tic
fX = fft(X,[],3);
%calculate the sum of the rank of block matrices from fft(bcirc(X0))
for i = 1:r
    Ti2 = fX(:,:,i);
    ras2(m) = ras2(m) + rank(Ti2);
end
t3(q,m) = toc
T3(m) = T3(m)+t3(q,m);

end
T1(m) = T1(m)/ite;
T3(m) = T3(m)/ite;
end
y = 1:le;
R = 2.^y;
loglog(R, T1,R,T3,'r-*')
% do the loglog plot
legend('prop1','prop3')
xlabel('r');
ylabel('time');

figure
plot(R, T1./R.^3,R,T3./R,'r-*')
legend('prop1','prop3')
xlabel('r');
ylabel('time');
T1./R.^3
T3./R
%to show T1 is O(r^3) and T3 is O(r)