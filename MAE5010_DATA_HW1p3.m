% Code by Christopher E. Petrin
% Prepared for MAE 5010: Data Assimilation (Spring 2020, Dr. Omer San)
% For Homework #1, Problem #3

% Generate a Cartesian grid: 0 <= x,y <= 1, dx = dy = 0.01
xgrid = 0:0.01:1;
ygrid = 0:0.01:1;

[X, Y] = meshgrid(xgrid, ygrid);
% Generate observation data at n random locations
%  Generate random locations
n = 2000;
zloc = rand(n,2);
figure(1)
plot(zloc(:,1), zloc(:,2),'+k')
xlim([0 1]);	ylim([0 1]);
%  Evaluate function at each pair
var = 0.1;
for i = 1:size(zloc,1)
	Z(i,1) = 2*zloc(i,1) + 4*zloc(i,2) + zloc(i,1)*zloc(i,2) + var*rand(1,1);
end

% Find i & j locations for each Z-point, so that the indices in the H
%  matrix may be known.
for c = 1:size(zloc,1)
	x = zloc(c,1);		
	xlow = floor(x*100); % Have to scale up by 100 in order to use floor()
	a = x - xlow/100; % Scale back down by 100
	y = zloc(c,2);		
	ylow = floor(y*100); % See above	
	b = y - ylow/100; % See above
	
	%i = find( xgrid*100 == xlow );
	% The above ( i = find() ) keeps failing for some reason... So reworked
	% simplistically below.
	i = xlow + 1; % Only works for given spacing...
	i_ = i + 1;
	
	%j = find( ygrid*100 == ylow );
	j = ylow + 1;
	j_ = j + 1;
	
	amat(c,1) = a; % each a coordinate is stored
	amat(c,2) = 0.01-a; % each abar coordinate is stored
	
	bmat(c,1) = b;
	bmat(c,2) = 0.01-b;
	
	% Column indices to locate stuff in the H matrix
	k(c,1) = (i-1)*size(xgrid,2)+j;
	k(c,2) = (i-1)*size(xgrid,2)+j_;
	k(c,3) = (i_-1)*size(xgrid,2)+j;
	k(c,4) = (i_-1)*size(xgrid,2)+j_;
end

% Populate/Compute H matrix
H = zeros(n,size(xgrid,2)*size(ygrid,2));
for c = 1:n
	H(c,k(c,1)) = amat(c,2)*bmat(c,2);
	H(c,k(c,2)) = amat(c,1)*bmat(c,2);
	H(c,k(c,3)) = amat(c,2)*bmat(c,1);
	H(c,k(c,4)) = amat(c,1)*bmat(c,1);
end

% Solve
xLU = f_xLU(H,Z);
xLUplot = reshape(xLU, 101,101);
figure(2)
contourf(X,Y,xLUplot, 'LineStyle','none')
colorbar
xlim([0 1]);	ylim([0 1]);

xQR = f_xQR(H,Z);
xQRplot = reshape(xQR, 101,101);
figure(3)
contourf(X,Y,xQRplot, 'LineStyle','none')
colorbar
xlim([0 1]);	ylim([0 1]);

xSVD = f_xSVD(H,Z);
xSVDplot = reshape(xSVD, 101,101);
figure(4)
contourf(X,Y,xSVDplot, 'LineStyle','none')
colorbar
xlim([0 1]);	ylim([0 1]);

normLU = norm(Z - H*xLU);
normQR = norm(Z - H*xQR);
normSVD = norm(Z - H*xSVD);