
files=codes/nrutil.c codes/nrutil.h codes/nr.h codes/mrqmin.c codes/mrqcof.c codes/covsrt.c codes/gasdev.c codes/gaussj.c codes/ran1.c codes/head.h codes/fit_poly.c codes/MC2.c codes/arrays.c

files+= codes/tcool.c codes/calculate_avg_values.c

main:
	gcc -o main main.c $(files) -lm

clean:
	@rm -f *.o main
	@echo "Object files cleaned."

