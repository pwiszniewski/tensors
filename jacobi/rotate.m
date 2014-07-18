function [v1 v2] = rotate (u1, u2, c, s)
	 
	 ## usage: [v1 v2] = rotate (u1, u2, c, s )
	 ## 
	 ## apply Jacobi rotation to a vector
	 ## v = J * u
	 ##
	 ## J = | c  s*  |
         ##     | -s  c* |
	 ##

 v1 = c*u1 + conj(s)*u2;
 v2 = -s*u1 + conj(c)*u2;

endfunction
