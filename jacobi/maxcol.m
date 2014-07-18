function l = maxcol ( a, k )
	 
	 ## usage: l = maxcol ( a, k )
	 ## 
	 ## find a maximal OFF-DIAGONAL element 
	 ## in a row k of matrix a to the right

  matsz = size(a,2);

  l = matsz;
  
  for j = k+1:matsz
      if ( abs(a(k,j)) > abs(a(k,l)) )
	 l = j;
      endif
  end

endfunction
