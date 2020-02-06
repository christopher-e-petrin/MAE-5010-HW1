% Code by Christopher E. Petrin
% Prepared for MAE 5010: Data Assimilation (Spring 2020, Dr. Omer San)
% For Homework #1, Problem #2

%% SET UP PROBLEM
% SETUP - Step 1: Pick four pairs (a_i, b_i) 1 <= i 1<= 4 of uniformly
%  distributed random numbers in range [0, 1]
a = rand(4,1);
b = rand(4,1);
%  Find pairs c and d from definition a + c = 1, b + d = 1
for i = 1:4
	c(i) = 1 - a(i);
	d(i) = 1 - b(i);
end

% SETUP - Step 2: Compute elements of the rows of H and verify that they
%  add up to 1.
%  Because problem statement in notes refers specifically to the matrix
%  generated in module 3.6, I am assuming that the indices of each value
%  remain the same (i.e. z1 is related to x1, x2, x5, and x6).
H1 = [c(1)*d(1) a(1)*d(1) 0 0 c(1)*b(1) a(1)*b(1) 0 0 0 0 0 0 0 0 0 0];
H2 = [0 0 c(2)*d(2) a(2)*d(2) 0 0 c(2)*b(2) a(2)*b(2) 0 0 0 0 0 0 0 0];
H3 = [0 0 0 0 0 c(3)*d(3) a(3)*d(3) 0 0 c(3)*b(3) a(3)*b(3) 0 0 0 0 0];
H4 = [0 0 0 0 0 0 0 0 0 0 c(4)*d(4) a(4)*d(4) 0 0 c(4)*b(4) a(4)*b(4)];

%  Create verification variables to be reported after; each should be equal
%  to 1.
sumH1 = sum(H1);		sumH2 = sum(H2);		sumH3 = sum(H3);

% SETUP - Step 3: Compute HH'
H = [H1; H2; H3; H4];
SPD = H*H';

% SETUP - Step 4: Generate observations Zi = 75 + Vi, Vi ~ N(0, sigma^2)
%  for 1 <= i <= 4; sigma^2 = 1.0 ( four observations)
vars = 1.0;
lowerbound = 0;		upperbound = vars;
for i = 1:4
	Z(i,1) = 75 + (lowerbound + (upperbound - lowerbound)*rand(1));
end

%% SOLVE PROBLEM
% PART A: LU
xLU = f_xLU(H,Z);
ZhLU = Z - H*xLU;
rxLU = Z - ZhLU;
norRLU = norm(ZhLU);

% PART B: QR
xQR = f_xQR(H,Z);
ZhQR = Z - H*xQR;
rxQR = Z - ZhQR;
norRQR = norm(ZhQR);

% PART C: SVD
xSVD = f_xSVD(H,Z);
ZhSVD = Z - H*xSVD;
rxSVD = Z - ZhSVD;
norRSVD = norm(ZhSVD);