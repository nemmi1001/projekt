CPP = g++
CPPFLAGS = -Wall -c

output: main.o wall.o library.o
	$(CPP) main.o wall.o library.o -o output

main.o: main.cpp
	$(CPP) $(CPPFLAGS) main.cpp

wall.o: wall.cpp wall.h
	$(CPP) $(CPPFLAGS) wall.cpp

knihovna.o: library.cpp library.h
	$(CPP) $(CPPFLAGS) library.cpp

clean:
	rm *.o  *.txt output


