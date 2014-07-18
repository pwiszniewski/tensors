function success = test_jacobi ()
	 
	 ## usage: success = test_jacobi ()
	 ## 
         ## apply different tests to test Jacobi  
	 ## diagonalization

  dim = 6;
  nmat = 4;

  A = zeros(dim,dim,nmat);
  B = zeros(dim,dim,nmat);

  switch T
    case 1
      ## test1 ## Diagonalize commuting hermitial matrices ##
    
      U = orth(rand(dim) + i*rand(dim));
      for j = 1:nmat
	A(:,:,j) = diag(rand(dim,1) + i*rand(dim,1));
	A(:,:,j) = U*A(:,:,j)*U';
      end
  
      [ V B ] = jacobi(A);
  
    case 2
      ## test2 ## Diagonalize slightly non-commuting Hermitean matrices ##

      U = orth(rand(dim) + i*rand(dim));
      sigma = 1e-3;
      
      for j = 1:nmat
	A(:,:,j) = diag(rand(dim,1) + i*rand(dim,1));
	K = rand(dim);
	K = 1/2 * (h + h');
	
	A(:,:,j) = expm(-i*sigma*K)*U*A(:,:,j)*U'*expm(i*sigma*K);
      end

      [ V B ] = jacobi(A);
  
  end
  
endfunction
