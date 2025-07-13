all: clean compile

clean:
	rm pqxt
	rm oqxt_eurosnp
	rm oqxt_eprint
	
compile:
	clear
	clear
	gcc -g -Wall pqxt.c -o pqxt -lcrypto -lm 
	gcc -g -Wall oqxt_eurosnp.c -o oqxt_eurosnp -lcrypto -lm 
	gcc -g -Wall oqxt_eprint.c -o oqxt_eprint -lcrypto -lm 

server:
	gdb ./server

client:
	gdb ./client
	
pqxt:
	gdb ./pqxt
	
oqxt_eurosnp:
	gdb ./oqxt_eurosnp
	
oqxt_eprint:
	gdb ./oqxt_eprint

run:
	./a.out

debug:
	gdb ./a.out
