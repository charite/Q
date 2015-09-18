
bin_PROGRAMS = $(top_builddir)/bin/Q

AM_CXXFLAGS=-lrt -I'seqan/core/include/' -Wall -Wno-long-long -pedantic -Wno-variadic-macros -fopenmp -static
AM_CXXFLAGS+=-DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 -O3 -DNDEBUG -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0
AM_LDFLAGS=-lbz2 -fopenmp
LDADD=-lz

__top_builddir__bin_Q_SOURCES = Q_EvaluateSignalsPvalue.cpp  Q_ReadInFiles.cpp Q.cpp Q_EvaluateSignalsPvalue.h \
Q_ReadInFiles.h Q_DetectSignals.cpp         Q_GetFragLength.cpp          Q_WriteResults.cpp \
Q_DetectSignals.h Q_GetFragLength.h Q_WriteResults.h Q_EvaluateSignalsCount.cpp  \
Q_ParseCommandLine.cpp Q_EvaluateSignalsCount.h Q_ParseCommandLine.h

EXTRA_DIST = seqan
