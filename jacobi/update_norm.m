function other_normf = update_norm (other_normf, diffB, k, l, key)
	 
	 ## usage: diff = update_norm (diffB, k, l, key)
	 ##        key - 0 - off-diagonal norm (default)
	 ##            - 1 - diagonal norm
	 ## 
	 ## Recalculates diagonal or off-diagonal Frobenius norm 
	 ## according to key
	 ##

  if (~exist('key','var'))
     key = 0;
  endif
	 
  matsz = size(diffB,2);
  
  if (key == 0)
    other_normf -= sum(vec(diffB));

    ## substract diagonal elements 
    other_normf += diffB(1,k) + diffB(2,l);
  else

    other_normf -= diffB(1,k) + diffB(2,l);
  endif

endfunction
