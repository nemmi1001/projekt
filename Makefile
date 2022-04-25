CPP = g++
CPPFLAGS = -Wall -c

output: main.o wall.o
	$(CPP) main.o wall.o -o output

main.o: main.cpp
	$(CPP) $(CPPFLAGS) main.cpp

wall.o: wall.cpp wall.h
	$(CPP) $(CPPFLAGS) wall.cpp

clean:
	rm *.o  *.txt output


