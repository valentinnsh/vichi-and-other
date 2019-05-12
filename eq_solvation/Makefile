files=$(wildcard *.f95)
obj=$(patsubst %.f95, %.o, $(files))

res: result 

result: main
main: $(obj)                                                          
	gfortran $^ -o $@
%.o: %.f95
	gfortran -c $< -o $@
clean:
	rm -f *.o main
