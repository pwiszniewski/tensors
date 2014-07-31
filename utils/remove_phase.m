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
	m(:,j) *= conj( m(2,j) / abs( m(2,j) ) );
    end
  elseif (mode == 'r')
    for j=1:nrow
      m(j,:) *= conj( m(j,2) / abs( m(j,2) ) );
    end
  else
      error('remove phase','wrong mode: %c',mode);
  endif
  
  mm = m;

endfunction
