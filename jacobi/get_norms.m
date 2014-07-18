function [normf other_normf] = get_norms (B, key)
	 
	 ## usage: [normf other_normf] = get_norms (B, key)
	 ##        key - 0 - off-diagonal norm (default)
	 ##            - 1 - diagonal norm
	 ##
	 ## computes total and off-diagonal (or diagonal)
         ## square of Frobenius norm of a matrix B. 
         ## The type of the computation is set by key.

  if (~exist('key','var'))
     key = 0;
  endif
	 
  normf = sum(vec(B));
  other_normf = sum(diag(B));
		     
  if (key == 0)
    other_normf = normf - other_normf;
  endif
  
endfunction
