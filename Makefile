all:
	gcc main.c -o convective_flow -lm
vg: all
	valgrind ./convective_flow
