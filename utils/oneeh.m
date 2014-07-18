function h = oneeh (dim,t,dopbc)
	 
	 ## usage: h = oneeh (dim)
	 ## 
	 ## initialize hopping matrix 
	 ## for Hubbard Hamiltonian
	 
  h = zeros(dim);
  
  if ( dopbc == 1 )
    for j = 1 : dim
      if ( j == 1 )
	h(j,j+1) = t;
	h(j,dim) = t;
	
      elseif ( j == dim )
	h(j,j-1) = t;
	h(j,1) = t;
	
      else
	h(j,j-1) = t;
	h(j,j+1) = t;
      end
    end

  else
    for j = 1 : dim
      if ( j == 1 )
	h(j,j+1) = t;
	
      elseif ( j == dim )
	h(j,j-1) = t;
	
      else
	h(j,j-1) = t;
	h(j,j+1) = t;
      end
    end
  end

endfunction
