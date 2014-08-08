function [ normd normof ] = bidiagnorm (EI, isym)
	 
	 ## usage: [ normd normof ] = bidiagnorm (EI)
	 ## 
	 ## 
	 
  dim = size(EI,1);
  
  normf = norm(vec(EI),'fro');

  normd = 0;

  if ( isym == 1 )

  ## Mulliken ordering

    for i = 1 : dim
      for j = 1 : dim
	if ( i != j )
	  normd += EI(i,i,j,j)*conj(EI(i,i,j,j));
	  normd += EI(i,j,j,i)*conj(EI(i,j,j,i));
	else
	  normd += EI(i,i,j,j)*conj(EI(i,i,j,j));
	endif
      end
    end
    
  elseif ( isym == 2 )
 
  ## Dirak ordering
      
    for i = 1 : dim
      for j = 1 : dim
	if ( i != j )
	  normd += EI(i,j,i,j)*conj(EI(i,j,i,j));
	  normd += EI(i,j,j,i)*conj(EI(i,j,j,i));
	else
	  normd += EI(i,j,i,j)*conj(EI(i,j,i,j));
	endif
      end
    end

  endif

  normof = normf^2 - normd;

endfunction
