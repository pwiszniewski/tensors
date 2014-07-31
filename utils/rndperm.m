function u = rndperm (dim)
	 
	 ## usage: u = rndperm (dim)
	 ## 
	 ## Generates a random permutation matrix with 
	 ## complex elements

	 a = rand(dim,1);
	 a = arrayfun(@exp, i*a);
	 perm = randperm(dim);
	 u = diag(a);
	 u = u(:,perm);

endfunction
