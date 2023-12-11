all:
	gcc main.c -o convective_flow -lm -Wall
vg: all
	valgrind ./convective_flow
