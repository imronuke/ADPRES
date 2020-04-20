module control

  use sdata, only: dp
  implicit none
  save

contains

  !******************************************************************************!

  SUBROUTINE GetPinPow()

    !
    ! Purpose:
    !    To perform pin power reconstruction
    !

    use sdata, only: nnod, f0, aprad, apaxi, afrad, ftem, mtem, cden, bcon, bpos
    use io,    only: AsmPow, AxiPow, AsmFlux, inp_read
    use xsec,  only: XS_updt
    use cmfd,  only: outer, powdis

    IMPLICIT NONE




  END SUBROUTINE GetPinPow




end module
