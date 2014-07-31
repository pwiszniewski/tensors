function v = triu2vec (t,diag)
	 
	 ## usage: v = triu2vec (t, diag)
	 ## 
	 ## upper triangular matrix to vector 
	 
  dim = size(t,1);
  patt = triu(ones(dim), diag);
  v = t(find(patt));

endfunction
