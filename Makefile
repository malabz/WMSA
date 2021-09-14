CC = gcc
#CC = icc
CFLAGS = -O3
INSTALLFOLDER = /usr/bin
CDHITFOLDER = $(INSTALLFOLDER)/cd-hit
MAFFTFOLDER = $(INSTALLFOLDER)/mafft
# You can edit this folder in order to install the program to another place

# Comment out the above line if your compiler 
# does not support TLS (thread-local strage).
ENABLE_MULTITHREAD = -Denablemultithread
THREADS = 1

ifdef ENABLE_MULTITHREAD
LIBS = -lm  -lpthread
else
LIBS = -lm
endif

ifdef ENABLE_ATOMIC
STDF = -std=c11
else
STDF = -std=c99
endif

MYCFLAGS = $(MNO_CYGWIN) $(ENABLE_MULTITHREAD) $(ENABLE_ATOMIC) $(STDF) $(CFLAGS)
HEADER = io.h msa.h

INSTALL = install

MSAMAINOBJ = mtxutl.o io.o constants.o msa_main.o function.o

wmsa : $(MSAMAINOBJ) 
	$(CC) -o $@ $(MSAMAINOBJ) $(MYCFLAGS) $(LIBS)

msa_main.o : msa_main.c $(HEADER)
	$(CC) $(MYCFLAGS) -c msa_main.c 

io.o : io.c $(HEADER)
	$(CC) $(MYCFLAGS) -c io.c 

mtxutl.o : ./mafft/mtxutl.c 
	$(CC) $(MYCFLAGS) -c ./mafft/mtxutl.c 

constants.o : constants.c
	$(CC) $(MYCFLAGS) -c constants.c

function.o : function.c function.h
	$(CC) $(MYCFLAGS) -c function.c

clean:
	rm *.o cd-hit/*.o mafft/*.o tmp* wmsa cd-hit/cd-hit cd-hit/cd-hit-454 cd-hit/cd-hit-est mafft/staralign mafft/disttbfast
	rm -rf swap

all:
	cd mafft && make -j$(THREADS)
	cd cd-hit && make -j$(THREADS)
	make -j$(THREADS)

install: 
	mkdir -p $(CDHITFOLDER) $(MAFFTFOLDER)
	cp wmsa $(INSTALLFOLDER)
	cp ./cd-hit/cd-hit ./cd-hit/cd-hit-est $(CDHITFOLDER)
	cp ./mafft/profilealign ./mafft/staralign $(MAFFTFOLDER)

uninstall:
	rm $(INSTALLFOLDER)/wmsa
	rm -rf $(CDHITFOLDER) $(MAFFTFOLDER)