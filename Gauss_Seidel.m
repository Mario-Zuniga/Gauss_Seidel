% Matrix X is generated
[x,y] = meshgrid(0 : pi / 250 : pi, 0 :pi / 250:pi);
X = cos(x.^2 + y.^2);

% Matrix H is generated
[f,c] = size(X);
H = randn(f,f);

% Matrix H is converted into a dominant
Temp = abs(H);
Temp2 = sum(Temp,2);
H = (H - diag(diag(H)) + diag(Temp2));

% Matrix X is distored
Y = H * X;

% A nex Matrix is generated for the filtered results of
% Gauss Seidel
XC = zeros(f,f);

% Tolerance is setled
tol = 1;

% In this cycle, the Gauss Seidel function is called
for i = 1 : f
   XC(:,i) = Gauss_Seidel_Funct(H, Y(:,i), tol); 
end

% Recovered image is shown
mesh(XC)


% Gauss Seidel function
function[xn] = Gauss_Seidel_Funct(A, b, tol)

% We need the number of columns & rows
[f,c] = size(A);

% Diagonal is obtained for the cycle
B = A - diag(diag(A));

error = 100;

xa = zeros(f,1);

% This cycle will continue as long as the error is
% greater than the tolerance
while error > tol
    
    % Previous value is saved
    xn = xa;
    
    % We obtain the values in the following cycles
    for i = 1 : f
        xn(i) = b (i);
        
        for j = 1 : c
           xn(i) = xn(i) - B(i,j) * xn(j);
        end
        
        xn(i) = xn(i) / A(i,i);
    end
    
    % Formula for the error
    error = norm(xn-xa) / norm(xn) * 100;
    xa = xn;
end
end