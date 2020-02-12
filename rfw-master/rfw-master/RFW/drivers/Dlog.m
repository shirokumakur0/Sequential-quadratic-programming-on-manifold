function L = Dlog(A,E)
% log_derivative: Derivative of the matrix logarithm function
%   log_derivative(A,E) is the derivative of the matrix logarithmic 
%   function at the point A in the direction of tangent vector E, 
%   calculated as described in [1].
%
%    [1] N.J. Higham, Functions of Matrices: Theory and Computation, 2008

n=size(A,1);
alpha=[3.1416*sqrt(-1),-3.1416*sqrt(-1),6.29*sqrt(-1),-6.29*sqrt(-1),1.0281e1*sqrt(-1),-1.0281e1*sqrt(-1),2.8894e1*sqrt(-1),-2.8894e1*sqrt(-1)];
beta=[1.5708*sqrt(-1),-1.5708*sqrt(-1),4.7125*sqrt(-1),-4.7125*sqrt(-1),7.9752*sqrt(-1),-7.9752*sqrt(-1),1.4823e1*sqrt(-1),-1.4823e1*sqrt(-1)];

s=0;
B=A;

while (norm(eye(n)-B,1)>(1-1/exp(1)))
    s=s+1;
    B=sqrtm(B);
end

E0=2^s*E;

for i=1:s
    E0=lyap(A^(0.5^i),A^(0.5^i),-E0);
end

G=lyap(B,B,-E0);

X=logm(B);

for i=8:-1:1
    G=lyap(eye(n)+X/alpha(i),eye(n)-X/alpha(i),-(eye(n)+X/(beta(i)))*G-G*(eye(n)-X/beta(i)));
end

L=2*real(G);




