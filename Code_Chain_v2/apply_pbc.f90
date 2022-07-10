!! This subroutines sets the periodic boundary conditions on
!! the padding elements
subroutine apply_pbc
  use mod_2dflu
  implicit none
  !% -- Apply PBC
  !! %-- ukx
  !% -- North and south faces
  ukx(1:n1, 0   ) = ukx(1:n1, n2);
  ukx(1:n1, n2+1) = ukx(1:n1,  1);
  ukx(1:n1, -1  ) = ukx(1:n1, n2-1);
  ukx(1:n1, n2+2) = ukx(1:n1,  2);


  !% --  East and west faces
  ukx(n1+1, 0:n2+1) = ukx(1 , 0:n2+1); 
  ukx(0   , 0:n2+1) = ukx(n1, 0:n2+1);
  ukx(n1+2, 0:n2+1) = ukx(2 , 0:n2+1); 
  ukx(-1   , 0:n2+1) = ukx(n1-1, 0:n2+1);

  !! %-- uky
  !% -- North and south faces
  uky(1:n1, 0   ) = uky(1:n1, n2);
  uky(1:n1, n2+1) = uky(1:n1,  1);
  uky(1:n1, -1  ) = uky(1:n1, n2-1);
  uky(1:n1, n2+2) = uky(1:n1,  2);


  !% --  East and west faces
  uky(n1+1, 0:n2+1) = uky(1 , 0:n2+1); 
  uky(0   , 0:n2+1) = uky(n1, 0:n2+1);
  uky(n1+2, 0:n2+1) = uky(2 , 0:n2+1); 
  uky(-1   , 0:n2+1) = uky(n1-1, 0:n2+1);
  !
end subroutine apply_pbc
