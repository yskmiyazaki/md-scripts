# definition of macro
CXX			= g++
CXXFLAGS	= -O3
LDFLAGS		=
INCLUDES	= -I.
TARGET		= bend_mod
OBJS		= main.o analyze.o coordinate.o fftsg.o

# generation rule
$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	-rm -f $(TARGET) $(OBJS) .nfs* *~ \#* core

.c.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<


main.o: analyze.h
analyze.o: analyze.h complex.h

