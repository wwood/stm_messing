all:
	gcc -g -Wall -O2 8mer_count.c -Iext -o 8mer_count -fgnu-tm -lz -lpthread
