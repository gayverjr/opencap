
#ifndef NUMGRID_EXPORT_H
#define NUMGRID_EXPORT_H

#ifdef NUMGRID_STATIC_DEFINE
#  define NUMGRID_EXPORT
#  define NUMGRID_NO_EXPORT
#else
#  ifndef NUMGRID_EXPORT
#    ifdef numgrid_shared_EXPORTS
        /* We are building this library */
#      define NUMGRID_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define NUMGRID_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef NUMGRID_NO_EXPORT
#    define NUMGRID_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef NUMGRID_DEPRECATED
#  define NUMGRID_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef NUMGRID_DEPRECATED_EXPORT
#  define NUMGRID_DEPRECATED_EXPORT NUMGRID_EXPORT NUMGRID_DEPRECATED
#endif

#ifndef NUMGRID_DEPRECATED_NO_EXPORT
#  define NUMGRID_DEPRECATED_NO_EXPORT NUMGRID_NO_EXPORT NUMGRID_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef NUMGRID_NO_DEPRECATED
#    define NUMGRID_NO_DEPRECATED
#  endif
#endif

#endif /* NUMGRID_EXPORT_H */
