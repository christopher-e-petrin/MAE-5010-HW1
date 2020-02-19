% Code by Christopher E. Petrin
% Prepared for MAE 5010: Data Assimilation (Spring 2020, Dr. Omer San)
% For Homework #1, Problem #4

%INPUT
clear
convPower = -5; % exponent used for 1*10^convPower to test for convergence

%Part A: Distribute two observations in each of 4x4 grid boxes for m = 18
% observations
ctr = 1;
for i = 1:4
	for j = 1:4
		xyloc(ctr,1) = i;		xyloc(ctr,2) = j;
		ctr = ctr + 1;
	end
end

zloc = zeros(1,2);
ctr = 1;
for i = 1:3
	for j = 1:3
		xys = rand(4,1);
		zloc(ctr,1) = xys(1) + i;
		zloc(ctr+1,1) = xys(2) + i;
		zloc(ctr,2) = xys(3) + j;
		zloc(ctr+1,2) = xys(4) + j;
		ctr = ctr + 2;
	end
end

% Plot grid & points to visually check for 2 points per box
figure(1)
plot(zloc(:,1),zloc(:,2),'ob') 
axis([1 4 1 4]);	xticks([1:4]); yticks([1:4]);	grid on

%PART B: Build the interpolation matrix H, (H is 18x16)
for ctr = 1:size(zloc,1)
	xn = zloc(ctr,1);		yn = zloc(ctr,2);
	
	a(ctr,1) = xn - floor(xn);
	a(ctr,2) = 1 - a(ctr,1);
	
	b(ctr,1) = yn - floor(yn);
	b(ctr,2) = 1 - b(ctr,1);
	
	i = floor(xn);
	i_ = i + 1;
	j = floor(yn);
	j_ = j + 1;
	
	% Find index values for a, abar, b, bbar
	k(ctr,1) = find( (xyloc(:,1) == i) & (xyloc(:,2) == j) );
	k(ctr,2) = find( (xyloc(:,1) == i) & (xyloc(:,2) == j_) );
	k(ctr,3) = find( (xyloc(:,1) == i_) & (xyloc(:,2) == j) );
	k(ctr,4) = find( (xyloc(:,1) == i_) & (xyloc(:,2) == j_) );
end
% Populate H matrix
H = zeros(18,16);
for ctr = 1:18
	H(ctr,k(ctr,1)) = a(ctr,2)*b(ctr,2);
	H(ctr,k(ctr,2)) = a(ctr,1)*b(ctr,2);
	H(ctr,k(ctr,3)) = a(ctr,2)*b(ctr,1);
	H(ctr,k(ctr,4)) = a(ctr,1)*b(ctr,1);
end

%PART C: Let Z = (z1, z2, ..., z18)' be the vector where zi = 70 + vi,
% where vi ~ N(0,variance)
var = 1;	% variance not specified in problem statement, setting to 1.
for i = 1:18
	Z(i,1) = 70 + rand(1);
end

%PART D: Construct f(x) = 0.5*[ x'(H'H)x - 2Z'Hx + Z'Z ]
% Guess x0 as the first row of H as calculated already. NO particular
% reason for doing this, just need a first guess.
xg(:,1) = H(1,:)';
xc(:,1) = H(1,:)';

% Construct parts of f(x) that aren't dependent on x.
%  f(x) = x'Ax - B'x + C
A = 0.5 * (H'*H);
B = Z'*H;
C = 0.5*Z'*Z;

%PART E: Apply Gradient & Conjugate gradient algorithm to minimize f(x)
rg(:,1) = B' - A*xg(:,1);
rc(:,1) = B' - A*xg(:,1);

% Conjugate algorithm
tic;
p(:,1) = rc(:,1);
for kc = 1:size(xc,1)
	%Step 0: Set up shorthad variables
	p_k = p(:,kc);	x_k = xc(:,kc);	r_k = rc(:,kc);
	
	%Step 1: Calculate alpha
	alpha(kc) = (p_k'*r_k)/(p_k'*A*p_k);
		a_k = alpha(kc);	%Shorthand
		
	%Step 2: Iterate x
	xc(:,kc+1) = x_k + a_k*p_k;
		x_kp1 = xc(:,kc+1);	%Shorthand
		
	%Step 3: Calculate residual
	rc(:,kc+1) = r_k - a_k * A * p_k;
		r_kp1 = rc(:,kc+1);	%Shorthand
		
	%Step 4: Test for convergence
	if r_kp1'*r_kp1 <= 1*10^(convPower)
		break
	else
		beta(kc) = (r_kp1'*r_kp1)/(r_k'*r_k);
		p(:,kc+1) = r_kp1 + beta(kc)*p_k;
	end
	
	clear p_k x_k r_k a_k x_kp1 r_kp1
end
tc = toc

% Gradient algorithm
tic;
kg = 1;
convergeTest = 0;
while convergeTest == 0
	alpha(kg) = ( rg(:,kg)' * rg(:,kg) ) / ( rg(:,kg)' * A * rg(:,kg) );
	xg(:,kg+1) = xg(:,kg) + alpha(:,kg)*rg(:,kg);
	
	% Test for convergence
	eg(:,kg) = xg(:,kg) - pinv(A)*B';
	if abs(mean(eg(:,kg))) <= 1*10^(convPower)
		break
	else % Update r(kg)
		rg(:,kg+1) = rg(:,kg) - alpha(kg)*A*rg(:,kg);
		
		kg = kg+1;
	end	
end
tg = toc

%PART F: Plot f(xk) vs k for each method & comment
for k = 1:size(xc,2)
	fxc(k,1) = k;
	fxc(k,2) = xc(:,k)'*A*xc(:,k) - B*xc(:,k) + C;
end

for k = 1:size(xg,2)
	fxg(k,1) = k;
	fxg(k,2) = xg(:,k)'*A*xg(:,k) - B*xg(:,k) + C;
end

for k = size(xc,2)+1:size(fxg,1)
	fxc(k,1) = k;
	fxc(k,2) = fxc(size(xc,2),2);
end

figure(2)
plot(fxc(1:100,1), fxc(1:100,2), '-b', 'LineWidth',2, 'DisplayName', 'Conjugate')
hold on
	plot(fxg(1:100,1), fxg(1:100,2), '-r', 'LineWidth',1, 'DisplayName', 'Gradient')
	legend('Location','Southeast')
	xlabel('k')
	ylabel('f(x_k)')
hold off
