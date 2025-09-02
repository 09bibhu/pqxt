all: clean compile

clean:
	rm pqxt
	rm oqxt_eurosnp
	
compile:
	clear
	clear
	gcc -g -Wall pqxt.c -o pqxt -lcrypto -lm 
	gcc -g -Wall oqxt_eurosnp.c -o oqxt_eurosnp -lcrypto -lm 

pqxt:
	gdb ./pqxt
	
oqxt_eurosnp:
	gdb ./oqxt_eurosnp
	
run:
	./a.out

debug:
	gdb ./a.out
