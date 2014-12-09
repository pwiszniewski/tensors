function EI = symmetrize (ET, isym, cplx)
	 
	 ## usage: EI = symmetrize (ET, isym, cplx)
	 ## symmetrizes integrals in the Mulliken or Dirak ordering 
	 ## 
         ##    cplx = 1 apply complex 4 fold symmetry
	 ##         = 0 apply real 8 fold symmetry
	 ##
         ##    isym = 0 - supersymm
	 ##           1 Mulliken
	 ##           2 Dirak 
	 ##
	 
  dim = size(ET,1);
  EI = zeros(dim,dim,dim,dim);
  
  if ( isym == 0 )
  ## cumulant symmetries [ij|kl] = [ji|lk]^* = [lk|ji] = [kl|ij]^* 
  ## real and complex cases are the same
     if ( cplx == 1 )
       for i = 1 : dim
	 for j = 1 : dim
	   for k = 1 : dim
	     for l = 1 : dim
	       EI(i,j,k,l) = (ET(i,j,k,l) + ET(l,j,k,i) + ET(l,k,j,i) + ET(i,k,j,l) 
			      + conj(ET(j,i,l,k)) + conj(ET(k,i,l,j)) + conj(ET(k,l,i,j)) + conj(ET(j,l,i,k))) / 8;		
	     end
	   end
	 end
       end
       
     else
       for i = 1 : dim
	 for j = 1 : dim
	   for k = 1 : dim
	     for l = 1 : dim
	       EI(i,j,k,l) = (ET(i,j,k,l) + ET(j,i,k,l) + ET(i,k,j,l)
			      + ET(j,k,i,l)+ ET(k,i,j,l) + ET(k,j,i,l)
			      + ET(i,j,l,k)+ ET(j,i,l,k) + ET(i,k,l,j)
			      + ET(j,k,l,i)+ ET(k,i,l,j) + ET(k,j,l,i)
			      + ET(i,l,j,k)+ ET(j,l,i,k) + ET(i,l,k,j)
			      + ET(j,l,k,i)+ ET(k,l,i,j) + ET(k,l,j,i)
			      + ET(l,i,j,k)+ ET(l,j,i,k) + ET(l,i,k,j)
			      + ET(l,j,k,i)+ ET(l,k,i,j) + ET(l,k,j,i)) / 24;		
	     end
	   end
	 end
       end
     endif
  elseif ( isym == 1 )
  ## symmetrize (ij|kl) = (kl|ij) = (ji|lk)^* = (lk|ji)^*

    if ( cplx == 1 )
      for i = 1 : dim
	for j = 1 : dim
	  for k = 1 : dim
	    for l = 1 : dim
	      EI(i,j,k,l) = (ET(i,j,k,l) + ET(k,l,i,j) + \
			     conj(ET(j,i,l,k)) + conj(ET(l,k,j,i)))/4;		
	    end
	  end
	end
      end
      
  ## symmetrize (ij|kl) = (kl|ij) = (ji|lk) = (lk|ji) =
  ## and also   (ji|kl) = (lk|ij) = (ij|lk) = (kl|ji) 

    else
      for i = 1 : dim
	for j = 1 : dim
	  for k = 1 : dim
	    for l = 1 : dim
	      EI(i,j,k,l) = (ET(i,j,k,l) + ET(j,i,k,l) 
			     + ET(i,j,l,k) + ET(j,i,l,k)
			     + ET(k,l,i,j) +  ET(l,k,i,j) 
			     + ET(k,l,i,j) + ET(l,k,j,i))/8;	
	    end
	  end
	end
      end
      
    endif
 
  elseif ( isym == 2 )

  ## symmetrize <ij|kl> = <ji|lk> = <kl|ij>^* = <lk|ji>^*

    if ( cplx == 1 )
      for i = 1 : dim
	for j = 1 : dim
	  for k = 1 : dim
	    for l = 1 : dim
	      EI(i,j,k,l) = (ET(i,j,k,l) + conj(ET(k,l,i,j)) + \
			     ET(j,i,l,k) + conj(ET(l,k,j,i)))/4;		
	    end
	  end
	end
      end
      
  ## symmetrize <ij|kl> = <ji|lk> = <kl|ij> = <lk|ji> =
  ## and also   <kj|il> = <li|jk> = <il|kj> = <jk|li> 

    else
      for i = 1 : dim
	for j = 1 : dim
	  for k = 1 : dim
	    for l = 1 : dim
	      EI(i,j,k,l) = (ET(i,j,k,l) + ET(j,i,l,k) + 
			     ET(i,l,k,j) + ET(l,i,j,k) + 
                             ET(k,j,i,l) + ET(j,k,l,i) +  
                             ET(k,l,i,j) + ET(l,k,j,i) ) / 8;

	    end
	  end
	end
      end
      
    endif

  endif

endfunction
