CHARM_PATH ?= ~/charm
CHARMC = $(CHARM_PATH)/bin/charmc

TARGET = fib

# P: number of PEs, N: fibonacci number
P = 2
N = 20

all: $(TARGET)

OBJS = $(TARGET).o

$(TARGET): $(OBJS)
	$(CHARMC) -language charm++ -o $@ $(OBJS)

$(TARGET).decl.h: $(TARGET).ci
	$(CHARMC) $<

$(TARGET).o: $(TARGET).C $(TARGET).decl.h
	$(CHARMC) -c $<

test: $(TARGET)
	./charmrun ++local ./$(TARGET) +p$(P) $(N)

clean:
	rm -f *.o *.decl.h *.def.h charmrun $(TARGET)
