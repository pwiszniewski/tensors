function t = vec2triu (v,dim,diag)
	 
	 ## usage: t = vec2triu (v)
	 ## 
	 ## packs a vector to an upper triangular matrix
	 
  t = triu(ones(dim),diag);
  ## add dimensionality checks here later
  t(t==1) = v;
  
endfunction
