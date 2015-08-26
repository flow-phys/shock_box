function setup_matrices()
global nA nB perx pery A1_L A1_U A2_L A2_U fil_alpha
global B1_L B1_U B2_L B2_U
%clear A1_U A1_L A2_U A2_L


%% A -> compact 6th order stuff
alpha = 1;  % The Diagonal weight
beta = 1/3;    % The off-diagonal weight

m = nA;
if (m>3)
A = diag( alpha*ones(m,1) )+ ...
    diag( beta*ones(m-1,1),1)+ ...
    diag( beta*ones(m-1,1),-1);
if (perx ==1) 
    A = A + diag( beta*ones(1,1),m-1);
    A = A + diag( beta*ones(1,1),-m+1);
else
    A(1,2) = 0;
    A(2,1) = 0;
    A(3,2) = 0;
    A(2,3) = 0;
    A(m,m-1) = 0;
    A(m-1,m) = 0;
    A(m-2,m-1) = 0;
    A(m-1,m-2) = 0;
end
[A1_L,A1_U] = lu(A);
end

clear A;
m = nB;
if (m>3)
A = diag( alpha*ones(m,1) )+ ... 
    diag( beta*ones(m-1,1),1)+ ...
    diag( beta*ones(m-1,1),-1);
if (pery ==1) 
    A = A + diag( beta*ones(1,1),m-1);
    A = A + diag( beta*ones(1,1),-m+1);
else
    A(1,2) = 0;
    A(2,1) = 0;
    A(3,2) = 0;
    A(2,3) = 0;
    A(m,m-1) = 0;
    A(m-1,m) = 0;
    A(m-2,m-1) = 0;
    A(m-1,m-2) = 0;
end
[A2_L,A2_U] = lu(A);
end



%% B -> compact 8th order filter
alpha = 1;  % The Diagonal weight
beta = fil_alpha;    % The off-diagonal weight

clear A;
m = nA;
if (m>5)
A = diag( alpha*ones(m,1) )+ ...
    diag( beta*ones(m-1,1),1)+ ...
    diag( beta*ones(m-1,1),-1);
if (perx ==1) 
    A = A + diag( beta*ones(1,1),m-1);
    A = A + diag( beta*ones(1,1),-m+1);
else
    A(1,2) = 0;
    A(2,1) = 0;
    A(2,3) = 0;
    A(3,2) = 0;
    A(3,4) = 0;
    A(4,3) = 0;
    A(4,5) = 0;
    A(m,m-1) = 0;
    A(m-1,m) = 0;
    A(m-1,m-2) = 0;
    A(m-2,m-1) = 0;
    A(m-2,m-3) = 0;
    A(m-3,m-2) = 0;
    A(m-3,m-4) = 0;
end
[B1_L,B1_U] = lu(A);
end

clear A;
m = nB;
if (m>3)
A = diag( alpha*ones(m,1) )+ ... 
    diag( beta*ones(m-1,1),1)+ ...
    diag( beta*ones(m-1,1),-1);
if (pery ==1) 
    A = A + diag( beta*ones(1,1),m-1);
    A = A + diag( beta*ones(1,1),-m+1);
else
    A(1,2) = 0;
    A(2,1) = 0;
    A(2,3) = 0;
    A(3,2) = 0;
    A(3,4) = 0;
    A(4,3) = 0;
    A(4,5) = 0;
    A(m,m-1) = 0;
    A(m-1,m) = 0;
    A(m-1,m-2) = 0;
    A(m-2,m-1) = 0;
    A(m-2,m-3) = 0;
    A(m-3,m-2) = 0;
    A(m-3,m-4) = 0;
end
[B2_L,B2_U] = lu(A);
end

end

