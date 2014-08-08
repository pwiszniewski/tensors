function [T w l] = init_sd_ten (dim)
	 
	 ## usage: T = init_sd_ten (dim)
	 ## 
	 ## initializes a trial 3-order tensor 
	 ## to test simultaneous diagonalization

  T = zeros(dim,dim,dim);

  w = rand(dim) + i*rand(dim);
  l = rand(dim) + i*rand(dim);

  for j = 1:dim
      T(:,:,j) = w * diag(l(j,:)) * w.';
  end

endfunction
