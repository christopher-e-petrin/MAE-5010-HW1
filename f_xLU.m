function [xLU] = f_xLU(H,Z)
	m = size(H,1);		n = size(H,2);
	if m < n
		xLU = H' * inv(H*H')*Z;
	elseif m > n
		xLU = inv(H'*H)*H'*Z;
	end
end

