main : main.o qwsim.o
	gcc  main.o qwsim.o -o main -lm

main.o : main.c qwsim.h
	gcc   -c main.c -lm

qwsim.o : qwsim.c qwsim.h
	gcc  -c qwsim.c -lm
