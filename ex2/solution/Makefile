CHARM_PATH ?= ~/charm
CHARMC = $(CHARM_PATH)/bin/charmc

TARGET = particle

# P: number of PEs, N: number of particles per cell, K: dimension of chare array
P = 4
N = 10000
K = 4

all: $(TARGET)

OBJS = $(TARGET).o

$(TARGET): $(OBJS)
	$(CHARMC) -O3 -language charm++ -o $@ $(OBJS)

$(TARGET).decl.h: $(TARGET).ci
	$(CHARMC) $<

$(TARGET).o: $(TARGET).C $(TARGET).decl.h
	$(CHARMC) -c $<

test: all
	./charmrun ++local ./$(TARGET) +p$(P) $(N) $(K)

clean:
	rm -f *.decl.h *.def.h conv-host *.o $(TARGET) charmrun *.sts *.gz *.projrc *.topo *.out
