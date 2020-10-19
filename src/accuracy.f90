module accuracy_m

    integer, parameter :: kint = 8
    integer, parameter :: kreal = kind(0.0d0)


    !big integer type digits
    integer, parameter :: desired_digits = 1000

    integer, parameter :: str_len_long = 2000
    real(kreal), parameter :: ccdist = 1.8d0

end module
