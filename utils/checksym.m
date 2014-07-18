function pass = checksym (EI, cplx, packing)
	 
	 ## usage: pass = checksym (EI, cplx, packing)
	 ## 
	 ## checks symmetries of the supplied EI
	 ## pass = 0 - all expected symmetries found
	 ## pass > 0 - some symmetries are broken. Bit wize
	 ##  0000001 - (ij|kl) = (kl|ij)   or <ij|kl> = <ji|lk>
         ##  0000010 - (ij|kl) = (ji|lk)^* or <ij|kl> = <kl|ij>^* 
         ##  0000100 - (ij|kl) = (lk|ji)^* or <ij|kl> = <lk|ji>^*
	 ## Additional real symmetries if cplx == 0: 
	 ##  0001000 - (ij|kl) = (ji|kl)   or <ij|kl> = <kj|il>
	 ##  0010000 - (ij|kl) = (lk|ij)   or <ij|kl> = <li|jk>
	 ##  0100000 - (ij|kl) = (ij|lk)   or <ij|kl> = <il|kj>
	 ##  1000000 - (ij|kl) = (kl|ji)   or <ij|kl> = <jk|li>

  dim = size(EI,1);
  
  pass = uint32(0);
  eps  = 1e-8; 

  if ( packing != 1 )

  ## check (ij|kl) = (kl|ij) = (ji|lk)^* = (lk|ji)^*

    if ( cplx == 1 )
      for i = 1 : dim
	for j = 1 : dim
	  for k = 1 : dim
	    for l = 1 : dim
	      if ( abs( EI(k,l,i,j) - EI(i,j,k,l) ) > eps )
		pass = bitset(pass,1);
	      endif
	      if ( abs( EI(j,i,l,k) - conj(EI(i,j,k,l)) ) > eps )
		 pass = bitset(pass,2);
	      endif
	      if ( abs( EI(l,k,j,i) - conj(EI(i,j,k,l)) ) > eps )
		 pass = bitset(pass,3);
	      endif
	    end
	  end
	end
      end
      
  ## check      (ij|kl) = (kl|ij) = (ji|lk) = (lk|ji) =
  ## and also   (ji|kl) = (lk|ij) = (ij|lk) = (kl|ij) 

    else

      for i = 1 : dim
	for j = 1 : dim
	  for k = 1 : dim
	    for l = 1 : dim
		if ( abs( EI(k,l,i,j) - EI(i,j,k,l) ) > eps)
		  pass = bitset(pass,1);
		endif
		if ( abs( EI(j,i,l,k) - EI(i,j,k,l) ) > eps)
		  pass = bitset(pass,2);
		endif
		if ( abs( EI(l,k,j,i) - EI(i,j,k,l) ) > eps )		
		  pass = bitset(pass,3);
		endif
		if ( abs( EI(j,i,k,l) - EI(i,j,k,l) ) > eps )
		  pass = bitset(pass,4);
		endif
		if ( abs( EI(l,k,i,j) - EI(i,j,k,l) ) > eps )		
		  pass = bitset(pass,5);
		endif
		if ( abs( EI(i,j,l,k) - EI(i,j,k,l) ) > eps )
		  pass = bitset(pass,6);
		endif
		if ( abs( EI(k,l,i,j) - EI(i,j,k,l) ) > eps )
		  pass = bitset(pass,7)
		endif
	    end
	  end
	end
      end
      
    endif

  else

  ## check <ij|kl> = <ji|lk> = <kl|ij>^* = <lk|ji>^*

    if ( cplx == 1 )
      for i = 1 : dim
	for j = 1 : dim
	  for k = 1 : dim
	    for l = 1 : dim
	      if ( abs( EI(j,i,l,k) - EI(i,j,k,l) ) > eps )
		pass = bitset(pass,1);
	      endif
	      if ( abs( EI(k,l,i,j) - conj(EI(i,j,k,l)) ) > eps )
		pass = bitset(pass,2);
	      endif
	      if ( abs( EI(l,k,j,i) - conj(EI(i,j,k,l)) ) > eps )
		 pass = bitset(pass,3);
	      endif
	    end
	  end
	end
      end
      
  ## check      <ij|kl> = <ji|lk> = <kl|ij> = <lk|ji> =
  ## and also   <kj|il> = <li|jk> = <il|kj> = <jk|li> 

    else

      for i = 1 : dim
	for j = 1 : dim
	  for k = 1 : dim
	    for l = 1 : dim
		if ( abs( EI(j,i,l,k) - EI(i,j,k,l) ) > eps )
		  pass = bitset(pass,1);
		endif
		if ( abs( EI(k,l,i,j) - EI(i,j,k,l) ) > eps )
		  pass = bitset(pass,2);
		endif
		if ( abs( EI(l,k,j,i) - EI(i,j,k,l) ) > eps )		
		  pass = bitset(pass,3);
		endif
		if ( abs( EI(k,j,i,l) - EI(i,j,k,l) ) > eps )
		  pass = bitset(pass,4);
		endif
		if ( abs( EI(l,i,j,k) - EI(i,j,k,l) ) > eps )		
		  pass = bitset(pass,5);
		endif
		if ( abs( EI(i,l,k,j) - EI(i,j,k,l) ) > eps )
		  pass = bitset(pass,6);
		endif
		if ( abs( EI(j,k,l,i) - EI(i,j,k,l) ) > eps )
		  pass = bitset(pass,7);
		endif
	    end
	  end
	end
      end
      
    endif

  endif
	 
endfunction
