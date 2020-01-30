function [xQR] = f_xQR(H,Z)
m = size(H,1);		n = size(H,2);
	if m < n
		[Q, R] = qr(H',0);
		xQR = Q*(inv(R')*Z);
	elseif m > n
		[Q, R] = qr(H,0);
		xQR = inv(R)*(Q'*Z);
	end
end

