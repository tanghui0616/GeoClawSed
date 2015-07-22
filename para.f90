
    !*****************************************************************
    !*
    !* This program  is used for set up precision
    !*
    !*****************************************************************
    module Set_Precision !Sets the default precision for floating point numbers

        implicit none

        public :: sngl, dbl, quad, Prec

        integer, parameter :: Sngl = Selected_Real_kind( 6, 37)
        integer, parameter :: Quad = Selected_Real_kind(33, 4931)
        integer, parameter :: Dbl  = Selected_Real_Kind(15, 307)!double precision
        integer, parameter :: Prec = Dbl ! Set default precision to double precision
        integer, parameter :: Charl = 30

    end module Set_Precision



















