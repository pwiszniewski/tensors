function t = multocum (t)
	 
	 ## usage: t = multocum (t)
	 ## 
	 ## reorder (special) mulliken tensor to cumulant
	 ## note that issues are possible if some complex symmetry is broken!!!!!!!!!!

  dim = size(t,1);

  for j = 1:dim
    for k = j+1:dim
	tmp = t(:,:,j,k);
	t(:,:,j,k) = t(:,:,k,j);
	t(:,:,k,j) = tmp;
    end
  end

endfunction
