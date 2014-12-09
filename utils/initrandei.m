function EI = initrandei (dim, rnk, isym, cplx, ee)
	 
	 ## usage: EI = initrandei (dim, rnk, isym, cplx, ee)
	 ## returns random integrals in the Mulliken or Dirak ordering 
	 ##  
	 ## isym = 0 supersymmetric cumulant
         ##        1 Mulliken
         ##        2 Dirak
	 ## ee - maximal norm difference

  if ( cplx == 1 )
    EI = rand(dim,dim,dim,dim) + i*rand(dim,dim,dim,dim);
  else
    EI = rand(dim,dim,dim,dim);
  endif

  ET = symmetrize(EI,isym,cplx);
  [EI trunc] = trunctensvd(ET,rnk);
  ET = symmetrize(EI,isym,cplx);

  count = 0;
  while ( (trunc > ee) && (count < 1000))

    [EI trunc]= trunctensvd(ET,rnk);
    ET = symmetrize(EI,isym,cplx);
    count += 1;

  endwhile
  
  printf("initrandei: %d iterations needed\n", count);
  
endfunction
