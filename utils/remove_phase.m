function mm = remove_phase (m,mode)
	 
	 ## usage: mm = remove_phase (m,mode)
	 ## 
	 ## Removes a phase from the rows or columns of 
         ## a matrix m as is determined by mode
	 ## 
	 ## m = 'c' - columns 
	 ##   = 'r' - rows

  [ nrow ncol ] = size(m);
  
  if (mode == 'c')
    for j=1:ncol

      for k=1:nrow
	if( abs(m(k,j)) > 1e-8 )
	  break
	end
      end

      m(:,j) *= conj( m(k,j) / abs( m(k,j) ) );
    end

  elseif (mode == 'r')
    for j=1:nrow
      
      for k=1:ncol
	if( abs(m(j,k)) > 1e-8 )
	  break
	end
      end
      
      m(j,:) *= conj( m(j,k) / abs( m(j,k) ) );
    end
  else
      error('remove phase','wrong mode: %c',mode);
  endif
  
  mm = m;

endfunction
