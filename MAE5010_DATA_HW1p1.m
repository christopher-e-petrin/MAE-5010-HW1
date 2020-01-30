% Code by Christopher E. Petrin
% Prepared for MAE 5010: Data Assimilation (Spring 2020, Dr. Omer San)
% For Homework #1, Problem #1

% Input Data, given in problem statement
f = [0.9; 1.0; 1.1; 1.2; 1.3];
gamf = [1/0.9; 1/0.7; 1/0.5; 1/0.3; 1/0.2];

% Discretized layers & temperatures, given in problem statement
p = [1.0; 0.5; 0.2];
T = [0.9 0.85 0.875];

% Set pre-noise temperatures in a variable
xbar = transpose(T);

% Evaluate ai1, ai2, ai3 using input data
%  Initialize outside of loop
a1 = zeros(5,1);		a2 = zeros(5,1);		a3 = zeros(5,1);

%  loop through all i indices, evaluating integral of function over
%   specified p.
for i = 1:5
	% Define function
	a_fun = @(p) p.*gamf(i).*exp(-p .* gamf(i));
	a1(i) = integral(a_fun, 0.5, 1.0);
	a2(i) = integral(a_fun, 0.2, 0.5);
	a3(i) = integral(a_fun, 0.0, 0.2);
end

% Collate into matrix H & compute general inverse of H
H = [a1, a2, a3];
Hplus = inv(H'*H)*H';

% Compute Zbar as Zbar = H*xbar
Zbar = H*xbar;

% Add noise to observation, then recover x
%  Given the following variances in problem statement
vars = [0.0; 0.1; 0.4; 0.8; 1.0; 1.2];
%  Initialize outside of loop
V = zeros(5,size(vars,1));
Z = zeros(5,size(vars,1));
xLU = zeros(3,size(vars,1),1);
xQR = zeros(3,size(vars,1),1);
LUresidual = zeros(size(vars,1),1);
QRresidual = zeros(size(vars,1),1);
SVDresidual = zeros(size(vars,1),1);

for j = 1:size(vars,1)
	V(:,j) = vars(j).*rand(5,1);
	% Let Z = Zbar + V be the noisy observation
	Z(:,j) = Zbar + V(:,j);
end

%% SOLVE THE PROBLEM
% Part A: Recover x using least squares method & plot results
tic
for j = 1:size(vars,1)
	% Recover x
	xLU(:,j) = Hplus*Z(:,j);
	%  Calculate residuals for each j variance index
	LUresidual(j) = norm(Z(:,j) - H*xLU(:,j));
end
LSresults = [vars LUresidual];
toc

% Part B: Recover x using QR decomposition
%  Use in-built MATLAB function for qr decomp. Note that, if using NVIDIA
%  drivers, this may crash MATLAB if drivers are out of date.
tic
[Q, R] = qr(H,0);
xQR = inv(R)*(Q'*Z);
%  Calculate residuals for each j variance index
for j = 1:size(vars,1)
		QRresidual(j) = norm(Z(:,j) - H*xQR(:,j));
end
QRresults = [vars QRresidual];
toc

% Part C: Recover x using SVD decomposition
tic
%  Step 1: Compute SVD of H
[U, S, V] = svd(H,'econ');
%  Step 2-4: Compute rotation, scaling & rotation
xSVD = V*(S^(-1))*U'*Z;
%  Calculate residuals for each j variance index
for j = 1:size(vars,1)
		SVDresidual(j) = norm(Z(:,j) - H*xQR(:,j));
end
SVDresults = [vars SVDresidual];
toc

% Plot
plot(LSresults(:,1), LSresults(:,2), 'ob','MarkerSize',10);
	xlabel('\sigma^2','FontName','TimesNewRoman','FontWeight','bold'); 
	ylabel('||Z - Hx_L_S||_2','FontWeight','bold')
	grid on
	hold on
	plot(QRresults(:,1), QRresults(:,2), 'xk', 'MarkerFaceColor','k',...
		'MarkerSize',10);
	plot(SVDresults(:,1), SVDresults(:,2), '+r', 'MarkerFaceColor','r',...
		'MarkerSize',10);
	hold off
