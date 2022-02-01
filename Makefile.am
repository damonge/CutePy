EXTRA_DIST = README.md

lib_LTLIBRARIES = libcute.la
libcute_la_SOURCES = src/cute.c
libcute_la_CPPFLAGS = $(OPENMP_CFLAGS) -I./src/
libcute_la_LDFLAGS = $(OPENMP_CFLAGS) -lcfitsio -lchealpix -lm

include_HEADERS = src/cute.h
