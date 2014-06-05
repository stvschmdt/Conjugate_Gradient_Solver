CC = gcc

CFLAGS = -Wall 

LDFLAGS = -lm

program = congrad

source = \
	 congrad.c \

obj = $(source:.c=.o)

$(program): $(obj) 
	        $(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.c
	        $(CC) $(CFLAGS) -c $< -o $@

clean:
	        rm -f $(program) $(obj)
run:
	        ./$(program)
edit:
	        vi -p $(source) 
