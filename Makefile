all:
	gcc -g -Wall main_parallel.c ext/threadpool.c -Iext -o stm_parallel -fgnu-tm -lz -lpthread
