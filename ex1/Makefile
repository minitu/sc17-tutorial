CHARM_PATH = ~/charm
CHARMC = $(CHARM_PATH)/bin/charmc

TARGET = fib

# P: number of PEs, N: fibonacci number
P = 4
N = 20

all: $(TARGET)

OBJS = $(TARGET).o

$(TARGET): $(OBJS)
	$(CHARMC) -language charm++ -o $@ $(OBJS)

$(TARGET).decl.h: $(TARGET).ci
	$(CHARMC) $<

$(TARGET).o: $(TARGET).C $(TARGET).decl.h
	$(CHARMC) -c $<

seq: $(TARGET)
	./charmrun ++local ./$(TARGET) +p1 $(N)

test: $(TARGET)
	./charmrun ++local ./$(TARGET) +p$(P) $(N)

clean:
	rm -f *.o *.decl.h *.def.h charmrun $(TARGET)
