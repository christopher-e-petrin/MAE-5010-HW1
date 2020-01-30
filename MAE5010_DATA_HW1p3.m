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
	x = zloc(c,1);		xlow = floor(x*100);	a = x - xlow/100;
	y = zloc(c,2);		ylow = floor(y*100);	b = y - ylow/100;
	
	%i = find( xgrid*100 == xlow );
	% The above keeps failing for some reason...
	i = xlow + 1; % Only works for given spacing...
	i_ = i + 1;
	
	%j = find( ygrid*100 == ylow );
	j = ylow + 1;
	j_ = j + 1;
	
	amat(c,1) = a; % each a coordinate is stored
	amat(c,2) = 0.01-a; % each abar coordinate is stored
	
	bmat(c,1) = b;
	bmat(c,2) = 0.01-b;
	
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
xLU = reshape(f_xLU(H,Z), 101,101);
figure(2)
contourf(X,Y,xLU, 'LineStyle','none')
colorbar
xlim([0 1]);	ylim([0 1]);

xQR = reshape(f_xQR(H,Z), 101,101);
figure(3)
contourf(X,Y,xQR, 'LineStyle','none')
colorbar
xlim([0 1]);	ylim([0 1]);

xSVD = reshape(f_xSVD(H,Z), 101,101);
figure(4)
contourf(X,Y,xSVD, 'LineStyle','none')
colorbar
xlim([0 1]);	ylim([0 1]);

norLU = norm(H*f_xLU(H,Z));
norQR = norm(H*f_xQR(H,Z));
norSVD = norm(H*f_xSVD(H,Z));
% for j = 1:size(Z,1)
% 	xLURes = norm(Z(j,1) - H*f_xLU(H,Z));
% 	xQRRes = norm(Z(j,1) - H*f_xQR(H,Z));
% 	xQRRes = norm(Z(j,1) - H*f_xSVD(H,Z));
% end