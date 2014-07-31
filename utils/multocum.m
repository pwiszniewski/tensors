function t = multocum (t)
	 
	 ## usage: t = multocum (t)
	 ## 
	 ## reorder mulliken tensor to cumulant

  dim = size(t,1);

  for j = 1:dim
    for k = j+1:dim
	tmp = t(:,:,j,k);
	t(:,:,j,k) = t(:,:,k,j);
	t(:,:,k,j) = tmp;
    end
  end

endfunction
