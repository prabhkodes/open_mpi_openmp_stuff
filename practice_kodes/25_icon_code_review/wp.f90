program wp_things
    implicit none

    real(kind=8) :: wp
    wp = 1.0d0

    print *, "the value of wp is->", wp
    print *, "kind(wp) is->", kind(wp)
    print *, "precision(wp) is->", precision(wp)
    print *, "range(wp) is->", range(wp)
    print *, "tiny(wp) is->", tiny(wp)
    print *, "huge(wp) is->", huge(wp)
    print *, "epsilon(wp) is->", epsilon(wp)
    print *, "radix(wp) is->", radix(wp)
    print *, "digits(wp) is->", digits(wp)
    print *, "maxexponent(wp) is->", maxexponent(wp)
    print *, "minexponent(wp) is->", minexponent(wp)
    print *, "storage_size(wp) is->", storage_size(wp)

end program wp_things